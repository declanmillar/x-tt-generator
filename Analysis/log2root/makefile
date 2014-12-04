# $@ name of target
# $^ name of all prerequisites with duplicates removed
# $< name of the first prerequisite

THEFILE = log2root

# Generate a list of C++ files that won't later be compiled to executables
SRC = .
LIB = .
EXECUTABLES = $(SRC)/$(THEFILE).cpp
CPPFILES = $(filter-out $(EXECUTABLES), $(wildcard $(SRC)/*.cpp))
# A list of object libraries, made by changing the file extensions in the above list
OBJS = $(patsubst $(SRC)/%.cpp, $(LIB)/%.o, $(CPPFILES))

# Header files
HEADERS = $(wildcard $(SRC)/*.h)

# Clang compiler w/ debug info, optimisation, warnings, and compile as a sharable option
CXX = clang++
# I would use -pedantic but it throws a billion warnings about the ROOT libraries
CXXFLAGS = -g -O -fPIC -Wall# -pedantic

# ROOT compiler flags
CXXFLAGS += $(shell root-config --cflags)

# Use Clang as the linker
LD = clang++

# ROOT linker flags, including the non-default TMVA, RooFit, RooStats, and NeuroBayes libraries
LDFLAGS = -O $(shell root-config --libs) -lTreePlayer

all: $(THEFILE)

$(THEFILE): $(LIB)/$(THEFILE).o
	@echo INFO: Linking $<
	$(LD) $^ $(LDFLAGS) -o $@

$(LIB)/%.o: $(SRC)/%.cpp
	@echo INFO: Compiling $<
	$(CXX) $(CXXFLAGS) -o $@ -c $<

.PHONY: clean cleanest tidy

clean:
	rm -f $(LIB)/*.o
	rm -f $(THEFILE)



# CXX             = clang
# # CXXFLAGS        = -O -Wall -fPIC -ggdb -Weffc++
# CXXFLAGS        = -O -Wall -fPIC -ggdb
# LD              = g++  
# ROOTCFLAGS      = -I$(ROOTSYS)/include  
# ROOTLIBS        = -L$(ROOTSYS)/lib  $(shell $(ROOTSYS)/bin/root-config --libs) 
# TREEPLAYER	:= -L $(ROOTSYS)/lib -lTreePlayer

# CXXFLAGS        += $(ROOTCFLAGS)   
# LIBS            = $(ROOTLIBS) $(TREEPLAYER)
# OBJS            =  log2root.o AtlasROOTStyle.o 

# PROGRAM         = log2root
  
# all:            $(PROGRAM)  
  
# $(PROGRAM):     $(OBJS)   
# 	                $(CXX) -g3 -pg -p -o $(PROGRAM) $(OBJS) $(LIBS) $(ROOTLIBS) $(CXXFLAGS) 
  
# clean:;         @rm -f $(OBJS) core *.o *.out   