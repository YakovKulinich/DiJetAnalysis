CXX = g++ 
CXXFLAGS = -Wall `root-config --cflags` -g -o2
LDFLAGS = -Wall `root-config --glibs` -Wl,-rpath=.

INCLUDEDIRS = -I./include

LIBDIRS     = -L.
LIBS        = -lMyLib -lTreePlayer

OBJ         = analysis
SOURCES     = DiJetAnalysis.cpp DiJetAnalysisData.cpp DiJetAnalysisMC.cpp \
	THmulf.cxx AtlasStyle.C MyRoot.C main.cpp
INCLUDES    = DiJetAnalysis.h

DSOURCES   = src/
DINCLUDES  = include/

MYLIBDIR    = lib/

CSOURCES    = $(addprefix $(DSOURCES),$(SOURCES))
CINCLUDES   = $(addprefix $(DINCLUDES),$(INCLUDES))

all: $(CSOURCES) $(CINCLUDES) libMyLib.so 
	$(CXX) $(CXXFLAGS) $(INCLUDEDIRS) $(CSOURCES) \
	$(LIBDIRS) $(LIBS) $(LDFLAGS) -o $(OBJ)

libMyLib.so: $(DINCLUDES)vTlv.h Dict.cpp
	g++ $(CXXFLAGS) Dict.cpp -shared -fPIC -o libMyLib.so 

Dict.cpp:
	rootcling -f Dict.cpp -rml libMyLib.so -c $(DINCLUDES)vTlv.h $(DINCLUDES)LinkDef.h

clean : 
	rm $(OBJ) Dict* libMyLib*
