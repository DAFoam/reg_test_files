#!/usr/bin/env bash

# the script will exit if there is any error
set -e

if [ -z "$WM_PROJECT" ]; then
  echo "OpenFOAM environment not found, forgot to source the OpenFOAM bashrc?"
  exit 1
fi

# for the remove command
shopt -s extglob 

# run pimple with SST
cp constant/turbulenceProperties_SST constant/turbulenceProperties
pimpleFoam

# rename U to UData
for t in $(seq 0.005 0.005 0.05)
do
   a=$(echo "$t" | sed 's/\0$//')
   cd $a && mv U.gz UData.gz && rm -rf -- !(UData.gz) && cd -
done

# decompose all times
decomposePar -time '0:'

# get ready for the SA run
cp constant/turbulenceProperties_SAFV3 constant/turbulenceProperties
rm -rf 0.*
