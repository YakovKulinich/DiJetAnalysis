CXX = g++ 
CXXFLAGS = -Wall `root-config --cflags` -g -o2
LDFLAGS = -Wall `root-config --glibs` -lboost_system -lboost_filesystem -Wl,-rpath=.

INCLUDEDIRS = -I./include

LIBDIRS     = -L.
LIBS        = -lMyLib -lTHmulf

OBJ         = analysis
SOURCES     = DiJetAnalysis.cpp DiJetAnalysisData.cpp DiJetAnalysisMC.cpp \
	      JetPair.cpp AtlasStyle.C MyRoot.cpp main.cpp
INCLUDES    = DiJetAnalysis.h

DSOURCES   = src/
DINCLUDES  = include/

MYLIBDIR    = lib/

CSOURCES    = $(addprefix $(DSOURCES),$(SOURCES))
CINCLUDES   = $(addprefix $(DINCLUDES),$(INCLUDES))

all: $(CSOURCES) $(CINCLUDES) libMyLib.so libTHmulf.so
	$(CXX) $(CXXFLAGS) $(CSOURCES) $(INCLUDEDIRS) \
	$(LIBDIRS) $(LIBS) $(LDFLAGS) -o $(OBJ)

libMyLib.so: $(DINCLUDES)dict.h Dict.cpp
	g++ $(CXXFLAGS) Dict.cpp -shared -fPIC -o libMyLib.so 

libTHmulf.so:$(DSOURCES)THmulf.cxx $(DINCLUDES)THmulf.h
	g++ `root-config --cflags` src/THmulf.cxx -I./include \
	-lTreePlayer `root-config --glibs` -shared -fPIC -o libTHmulf.so

Dict.cpp:
	rootcling -f Dict.cpp -rml libMyLib.so -c $(DINCLUDES)dict.h $(DINCLUDES)LinkDef.h

clean : 
	rm $(OBJ) Dict* libMyLib* libTHmulf*
