#!/bin/bash

# --evelibs needed for TGeoNode and TGeoManager
# Link directly to edep-sim
g++ -g -Ofast -o dumpSSRITree.exe dumpSSRITree.cpp \
  $(root-config --cflags --evelibs) \
  -I${EDEP_SIM}/include -L ${EDEP_SIM}/lib -ledepsim_io
