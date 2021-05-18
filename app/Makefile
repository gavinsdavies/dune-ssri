CXX=g++
AR=ar
LD_SHARED=g++

# Need evelibs for TGeoNode and TGeoManager
ROOT_LIBS = $(shell root-config --evelibs)
ROOT_INC = $(shell root-config --cflags)

EDEP_LIBS = -L $(EDEP_SIM)/lib -ledepsim_io
# Need CLHEP for unit conversion, can probably just move this into a standalone header to remove dependency?
EDEP_INC = -I$(EDEP_SIM)/include -I$(CLHEP_INC)

# The TMS includes
TMS_INC = -I../src

CXXFLAGS = -Wall -O3 -g -fPIC -Werror
CXXFLAGS += $(ROOT_INC) $(EDEP_INC) $(TMS_INC)

# Our lovely collection of dependencies
LIB_OBJ = $(EDEP_LIBS) $(ROOT_LIBS)

LIB_DIR = ../lib
LIB_TMS = -lTMS_Prod

all: ConvertToTMSTree.exe BetheBloch_Example.exe

ConvertToTMSTree.exe: ConvertToTMSTree.cpp
	$(CXX) $(CXXFLAGS) -o $@ ConvertToTMSTree.cpp -L $(LIB_DIR) $(LIB_TMS) $(LIB_OBJ) $(BUILD_FLAGS)

BetheBloch_Example.exe: BetheBloch_Example.cpp
	$(CXX) $(CXXFLAGS) $(ROOT_INC) -o $@ $^ $(ROOT_LIBS)

clean:
	rm -f lib* *.o *.a *.exe
