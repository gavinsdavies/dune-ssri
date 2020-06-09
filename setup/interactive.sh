#! /usr/bin/env bash

MY_SSRI_DIR=/dune/app/users/gsdavies/ssri

# edep-sim needs to know where a certain GEANT .cmake file is...
G4_cmake_file=`find ${GEANT4_FQ_DIR}/lib64 -name 'Geant4Config.cmake'`
export Geant4_DIR=`dirname $G4_cmake_file`

# edep-sim needs to have the GEANT bin directory in the path
export PATH=$PATH:$GEANT4_FQ_DIR/bin


export GXMLPATH=$MY_SSRI_DIR/gevgen:${GXMLPATH}
export GNUMIXML="GNuMIFlux.xml"

export LD_LIBRARY_PATH=$MY_SSRI_DIR/edep-sim/edep-gcc-6.4.0-x86_64-pc-linux-gnu/lib:${LD_LIBRARY_PATH}
export PATH=$MY_SSRI_DIR/edep-sim/edep-gcc-6.4.0-x86_64-pc-linux-gnu/bin:${PATH}
