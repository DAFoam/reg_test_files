#!/usr/bin/env python
"""
DAFoam run script for the NACA0012 airfoil at low-speed
"""

# =============================================================================
# Imports
# =============================================================================
import os
import argparse
from mpi4py import MPI
from dafoam import PYDAFOAM, optFuncs
from pygeo import *
from pyspline import *
from idwarp import USMesh
from pyoptsparse import Optimization, OPT
import numpy as np


# =============================================================================
# Input Parameters
# =============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("--opt", help="optimizer to use", type=str, default="slsqp")
parser.add_argument("--task", help="type of run to do", type=str, default="runAdjoint")
args = parser.parse_args()
gcomm = MPI.COMM_WORLD


TRef = 0.999

# Set the parameters for optimization
daOptions = {
    "designSurfaces": ["upperWall"],
    "solverName": "DAScalarTransportFoam",
    "adjPCLag": 1000,
    "adjJacobianOption": "JacobianFree",
    "printIntervalUnsteady": 1,
    "primalBC": {
        "T0": {"variable": "T", "patches": ["inlet"], "value": [TRef]}
    },
    "unsteadyAdjoint": {"mode": "timeAccurateAdjoint", "nTimeInstances": 3},
    "objFunc": {
        "TVOL": {
            "part1": {
                "type": "variableVolSum",
                "source": "boxToCell",
                "min": [-50.0, -50.0, -50.0],
                "max": [50.0, 50.0, 50.0],
                "varName": "T",
                "varType": "scalar",
                "component": 0,
                "isSquare": 0,
                "scale": 1.0,
                "addToAdjoint": True,
            }
        },
    },
    "debug": True,
    "primalMinResTol": 1e-16,
    #"adjStateOrdering": "cell",
    "adjEqnOption": {"pcFillLevel": 0, "jacMatReOrdering": "natural", "useNonZeroInitGuess": False},
    "normalizeStates": {"U": 1.0, "p": 1.0, "nuTilda": 0.1, "phi": 1.0},
    "adjPartDerivFDStep": {"State": 1e-7, "FFD": 1e-2},
    "designVar": {},
}

# mesh warping parameters, users need to manually specify the symmetry plane and their normals
meshOptions = {
    "gridFile": os.getcwd(),
    "fileType": "openfoam",
    # point and normal for the symmetry plane
    "symmetryPlanes": [[[0.0, 0.0, 0.0], [0.0, 0.0, -0.0005]], [[0.0, 0.0, 0.1], [0.0, 0.0, 0.0005]]],
}

# options for optimizers
if args.opt == "snopt":
    optOptions = {
        "Major feasibility tolerance": 1.0e-7,
        "Major optimality tolerance": 1.0e-7,
        "Function precision": 1.0e-7,
        "Verify level": -1,
        "Major iterations limit": 50,
        "Nonderivative linesearch": None,
        "Print file": "opt_SNOPT_print.out",
        "Summary file": "opt_SNOPT_summary.out",
    }
elif args.opt == "slsqp":
    optOptions = {
        "ACC": 1.0e-7,
        "MAXIT": 50,
        "IFILE": "opt_SLSQP.out",
    }
else:
    print("opt arg not valid!")
    exit(0)


# =============================================================================
# Design variable setup
# =============================================================================
DVGeo = DVGeometry("./FFD/bumpFFD.xyz")
# nTwists is the number of FFD points in the spanwise direction
nTwists = DVGeo.addRefAxis("bodyAxis", xFraction=0.25, alignIndex="k")
# select points
iVol = 0
pts = DVGeo.getLocalIndex(iVol)
indexList = pts[3, 0, 0].flatten()
PS = geo_utils.PointSelect("list", indexList)
# shape
#DVGeo.addGeoDVLocal("shapey", lower=-1.0, upper=1.0, axis="y", scale=1.0, pointSelect=PS)
#daOptions["designVar"]["shapey"] = {"designVarType": "FFD"}

def tin(val, geo):
    inletT = float(val[0])
    DASolver.setOption("primalBC", {"T0": {"variable": "T", "patches": ["inlet"], "value": [inletT]}})
    DASolver.updateDAOption()

DVGeo.addGeoDVGlobal("tbc", [TRef], tin, lower=0.0, upper=100.0, scale=1.0)
daOptions["designVar"]["tbc"] = {"designVarType": "BC", "patches": ["inlet"], "variable": "T", "comp": 0}
# =============================================================================
# DAFoam initialization
# =============================================================================
DASolver = PYDAFOAM(options=daOptions, comm=gcomm)
DASolver.setDVGeo(DVGeo)
mesh = USMesh(options=meshOptions, comm=gcomm)
DASolver.addFamilyGroup(DASolver.getOption("designSurfaceFamily"), DASolver.getOption("designSurfaces"))
DASolver.printFamilyList()
DASolver.setMesh(mesh)
evalFuncs = []
DASolver.setEvalFuncs(evalFuncs)

# =============================================================================
# Constraint setup
# =============================================================================
DVCon = DVConstraints()
DVCon.setDVGeo(DVGeo)
DVCon.setSurface(DASolver.getTriangulatedMeshSurface(groupName=DASolver.getOption("designSurfaceFamily")))

# =============================================================================
# Initialize optFuncs for optimization
# =============================================================================
def setObjFuncsUnsteady(DASolver, funcs, evalFuncs):
    nTimeInstances = DASolver.getOption("unsteadyAdjoint")["nTimeInstances"]
    for func in evalFuncs:
        avgObjVal = 0.0
        for i in range(1, nTimeInstances):
            avgObjVal += DASolver.getTimeInstanceObjFunc(i, func)
        funcs[func] = avgObjVal # / (nTimeInstances - 1)

    funcs["fail"] = False


def setObjFuncsSensUnsteady(DASolver, funcs, funcsSensAllTimeInstances, funcsSensCombined):

    nTimeInstances = 1.0 * len(funcsSensAllTimeInstances)
    for funcsSens in funcsSensAllTimeInstances:
        for objFunc in funcsSens:
            if objFunc != "fail":
                funcsSensCombined[objFunc] = {}
                for dv in funcsSens[objFunc]:
                    funcsSensCombined[objFunc][dv] = np.zeros_like(funcsSens[objFunc][dv], dtype="d")

    for funcsSens in funcsSensAllTimeInstances:
        for objFunc in funcsSens:
            if objFunc != "fail":
                for dv in funcsSens[objFunc]:
                    funcsSensCombined[objFunc][dv] += funcsSens[objFunc][dv] # / nTimeInstances

    funcsSensCombined["fail"] = False

    if gcomm.rank == 0:
        print(funcsSensCombined)
    return


optFuncs.DASolver = DASolver
optFuncs.DVGeo = DVGeo
optFuncs.DVCon = DVCon
optFuncs.evalFuncs = evalFuncs
optFuncs.gcomm = gcomm
optFuncs.setObjFuncsUnsteady = setObjFuncsUnsteady
optFuncs.setObjFuncsSensUnsteady = setObjFuncsSensUnsteady

# =============================================================================
# Task
# =============================================================================
if args.task == "opt":

    optProb = Optimization("opt", objFun=optFuncs.calcObjFuncValues, comm=gcomm)
    DVGeo.addVariablesPyOpt(optProb)
    DVCon.addConstraintsPyOpt(optProb)

    optProb.addObj("TVOL", scale=1)
    # optProb.addCon("CL", lower=CL_target, upper=CL_target, scale=1)

    if gcomm.rank == 0:
        print(optProb)

    DASolver.runColoring()

    opt = OPT(args.opt, options=optOptions)
    histFile = "./%s_hist.hst" % args.opt
    sol = opt(optProb, sens=optFuncs.calcObjFuncSens, storeHistory=histFile)
    if gcomm.rank == 0:
        print(sol)

elif args.task == "runPrimal":
#    xDV = DVGeo.getValues()
#    xDV["tbc"][0] = 0.01
#    optFuncs.calcObjFuncValuesUnsteady(xDV)
    optFuncs.runPrimal(objFun=optFuncs.calcObjFuncValuesUnsteady)

elif args.task == "runAdjoint":

    optFuncs.runAdjoint(objFun=optFuncs.calcObjFuncValuesUnsteady, sensFun=optFuncs.calcObjFuncSensUnsteady)

elif args.task == "verifySens":

    optFuncs.verifySens(objFun=optFuncs.calcObjFuncValuesUnsteady, sensFun=optFuncs.calcObjFuncSensUnsteady)

else:
    print("task arg not found!")
    exit(0)
