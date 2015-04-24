CXX             = g++  
# CXXFLAGS        = -O -Wall -fPIC -ggdb -Weffc++
CXXFLAGS        = -O -Wall -fPIC -ggdb
LD              = g++  
# ROOTCFLAGS      = -I$(ROOTSYS)/include  
ROOTCFLAGS      = $(shell root-config --cflags)
# ROOTLIBS        = -L$(ROOTSYS)/lib  $(shell $(ROOTSYS)/bin/root-config --libs) 
ROOTLIBS        = $(shell root-config --libs)
# TREEPLAYER	:= -L $(ROOTSYS)/lib -lTreePlayer

CXXFLAGS        += $(ROOTCFLAGS)   
LIBS            = $(ROOTLIBS) $(TREEPLAYER)
OBJS            =  main.o atlas_style.o RootTuple.o analysis_transverse.o 

PROGRAM         = go
  
all:            $(PROGRAM)  
  
$(PROGRAM):     $(OBJS)   
	                $(CXX) -g3 -pg -p -o $(PROGRAM) $(OBJS) $(LIBS) $(ROOTLIBS) $(CXXFLAGS) 
  
clean:;         @rm -f $(OBJS) core *.o *.out   
