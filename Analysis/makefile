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
OBJS            =  main.o atlas_style.o superpose.o superpose2d.o RootTuple.o analysis.o 

BIN = go
LIB = Library
SRC = Source
OUT = .

# Compile all files ending in .f in SRC
$(LIB)/%.o: $(SRC)/%.cpp 
	$(CXX) $(CXXFLAGS) -c -o  $@ $<

# Link mainfile and all processes
$(OUT)/$(BIN): $(patsubst %, $(LIB)/%, $(OBJS))
	$(LD) $(ROOTLIBS) -o $@ $^

# Clean up
clean:
	rm -f $(LIB)/*.o $(OUT)/$(BIN)
  
# all:            $(PROGRAM)  
  
# $(PROGRAM):     $(LIB)/$(OBJS)   
# 	                $(SRC)/$(CXX) -g3 -pg -p -o $(PROGRAM) $(LIB)/$(OBJS) $(LIBS) $(ROOTLIBS) $(CXXFLAGS) 
  
# clean:;         @rm -f $(OBJS) core *.o *.out   
