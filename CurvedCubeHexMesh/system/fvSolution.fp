/*--------------------------------*- C++ -*---------------------------------*\ 
| ========                 |                                                 | 
| \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | 
|  \    /   O peration     | Version:  v1812                                 | 
|   \  /    A nd           | Web:      www.OpenFOAM.com                      | 
|    \/     M anipulation  |                                                 | 
\*--------------------------------------------------------------------------*/ 
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

SIMPLE
{
    nNonOrthogonalCorrectors           0;
}

solvers
{
    p
    {
        
        solver                         GAMG;
        smoother                       GaussSeidel;
        relTol                         0.1;
        tolerance                      0;
        maxIter                        10;
    }
    "U|nuTilda"
    {
        solver                         smoothSolver;
        smoother                       GaussSeidel;
        relTol                         0.1;
        tolerance                      0;
        nSweeps                        1;
        maxIter                        10;
    }
    "pseudo_U|pseudoU"
    {
        solver                         smoothSolver;
        smoother                       GaussSeidel;
        relTol                         1e-2;
        tolerance                      0;
        nSweeps                        1;
        maxIter                        100;
    }
    "pseudo_p|pseudoP"
    {
        
        solver                         GAMG;
        smoother                       GaussSeidel;
        relTol                         1e-2;
        tolerance                      0;
        maxIter                        100;
    }
    "pseudo_nuTilda|pseudoNuTilda"
    {
        
        solver                         GAMG;
        smoother                       GaussSeidel;
        relTol                         1e-2;
        tolerance                      0;
        maxIter                        100;
    }
}

relaxationFactors
{
    fields
    {
        "(p|p_rgh|pseudo_p|pseudoP)"                    0.30;
    }
    equations
    {
        "(U|T|e|h|nuTilda|k|epsilon|omega|pseudo_U|pseudo_nuTilda|pseudoU|pseudoNuTilda)" 0.70;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
