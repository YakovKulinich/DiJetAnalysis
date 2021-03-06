ROOUNFOLDDIR = /home/yakov/FC/RooUnfold-1.1.1

CXX      = g++
CXXFLAGS = -Wall `root-config --cflags` -g -o2
LDFLAGS  = -Wall `root-config --glibs` -lboost_system -lboost_filesystem -Wl,-rpath=.

INCLUDEDIRS = -I./include

LIBDIRS     = -L.
LIBS        = -lMyLib

OBJ         = analysis
SOURCES     = DiJetAnalysis.cpp DiJetAnalysisData.cpp DiJetAnalysisMC.cpp \
	      DiJetAnalysisBoth.cpp DeltaPhiProj.cpp UncertaintyProvider.cpp \
	      HIJESUncertaintyProvider.cpp AtlasStyle.C  MyRoot.cpp main.cpp
INCLUDES    = DiJetAnalysis.h

DSOURCES    = src/
DINCLUDES   = include/

MYLIBDIR    = lib/

CSOURCES    = $(addprefix $(DSOURCES),$(SOURCES))
CINCLUDES   = $(addprefix $(DINCLUDES),$(INCLUDES))

all: $(CSOURCES) $(CINCLUDES) libMyLib.so 
	$(CXX) $(CXXFLAGS) $(CSOURCES) $(INCLUDEDIRS) \
	$(LIBDIRS) $(LIBS) $(LDFLAGS) -o $(OBJ)

libMyLib.so: $(DINCLUDES)dict.h Dict.cpp
	g++ $(CXXFLAGS) Dict.cpp -shared -fPIC -o libMyLib.so 

Dict.cpp:
	rootcling -f Dict.cpp -rml libMyLib.so -c $(DINCLUDES)dict.h $(DINCLUDES)LinkDef.h

clean : 
	rm $(OBJ) Dict* libMyLib* 
