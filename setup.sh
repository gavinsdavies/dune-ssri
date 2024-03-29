#!/bin/bash

# Example of setup needed for SSRI dependencies
# Here we use our own version of edep-sim, but you should be able to use one set up by ups

MY_SSRI_DIR=/dune/app/users/cwret/SSRI

# Copied from Gavin's SSRI setup
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup dk2nu        v01_05_01b   -q e15:prof
setup genie        v2_12_10c    -q e15:prof
setup genie_xsec   v2_12_10     -q DefaultPlusValenciaMEC
setup genie_phyopt v2_12_10     -q dkcharmtau
setup geant4       v4_10_3_p01b -q e15:prof
setup ifdhc

# edep-sim needs to find geant4.sh to build
export PATH=$PATH:${GEANT4_FQ_DIR}/bin

# edep-sim also needs cmake for setup
#setup cmake v3_9_0 -f Linux64bit+2.6-2.12
export GXMLPATH=$MY_SSRI_DIR/gevgen:${GXMLPATH}
export GNUMIXML="GNuMIFlux.xml"

# Finally Gavin's code needs
export EDEP_SIM=${MY_SSRI_DIR}/edep-sim/edep-gcc-6.4.0-x86_64-pc-linux-gnu
export LD_LIBRARY_PATH=${EDEP_SIM}/lib:${LD_LIBRARY_PATH}
export PATH=${EDEP_SIM}/bin:${PATH}

# Add TMS execs and library directory to env
export PATH=${PATH}:${PWD}/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/lib

echo "Setup TMS environment"
