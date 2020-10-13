CXX=g++
AR=ar
LD_SHARED=g++

# Need evelibs for TGeoNode and TGeoManager
ROOT_LIBS = $(shell root-config --evelibs)
ROOT_INC = $(shell root-config --cflags)

EDEP_LIBS = -L $(EDEP_SIM)/lib -ledepsim_io
EDEP_INC = -I$(EDEP_SIM)/include

CXXFLAGS = -Wall -Og -g -fPIC -Werror
CXXFLAGS += $(ROOT_INC) $(EDEP_INC)

# Replace and create
ARFLAGS = rc

# Get the ROOT CXX flags in the shared lib
SO_FLAGS = -shared $(ROOT_INC)

# Include this directory when looking for libs
BUILD_FLAGS = -Wl,-R -Wl,.

# Our lovely collection of libs
LIB_OBJ = $(EDEP_LIBS) $(ROOT_LIBS)

TMS_OBJ = TMS_Bar.o TMS_Hit.o TMS_TrueHit.o TMS_Event.o TMS_TrueParticle.o TMS_EventViewer.o

%.o:%.cpp
	$(CXX) -c $(CXXFLAGS) $<

all: libTMS_Prod ConvertToTMSTree.exe

libTMS_Prod: libTMS_Prod.a libTMS_Prod.so

ConvertToTMSTree.exe: ConvertToTMSTree.cpp libTMS_Prod
	$(CXX) $(CXXFLAGS) -o ConvertToTMSTree.exe ConvertToTMSTree.cpp -L . -lTMS_Prod $(LIB_OBJ) -Wl,-R -Wl,.

# Make the archive or static lib
libTMS_Prod.a: $(TMS_OBJ)
	$(AR) $(ARFLAGS) $@ $^

# Make the shared object
libTMS_Prod.so: libTMS_Prod.a $(TMS_OBJ)
	$(LD_SHARED) $(SO_FLAGS) $^ -o $@ $(LIB_OBJ)

clean:
	rm -f lib* *.o
