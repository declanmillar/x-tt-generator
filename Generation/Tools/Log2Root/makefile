# ROOT makefile
# Declan Millar <d.millar@soton.ac.uk>
BIN = log2root

# Generate a list of C++ files that won't later be compiled to executables
SRC = .
LIB = .
EXECUTABLES = $(SRC)/$(BIN).cpp
CPPFILES = $(filter-out $(EXECUTABLES), $(wildcard $(SRC)/*.cpp))
# A list of object libraries, made by changing the file extensions in the above list
OBJS = $(patsubst $(SRC)/%.cpp, $(LIB)/%.o, $(CPPFILES))

# Header files
HEADERS = $(wildcard $(SRC)/*.h)

# Clang compiler w/ debug info, optimisation, warnings, and compile as a sharable option
CPP = clang++
# I would use -pedantic but it throws a billion warnings about the ROOT libraries
CXXFLAGS = -g -O -fPIC -Wall #-pedantic -I/usr/local/Cellar/boost/1.57.0/include

# ROOT compiler flags
CXXFLAGS += $(shell root-config --cflags)

# Use Clang as the linker
LD = clang++

# ROOT linker flags, including the non-default TMVA, RooFit, RooStats, and NeuroBayes libraries
LDFLAGS = -O $(shell root-config --libs) -lTreePlayer

all: $(BIN)

$(BIN): $(LIB)/$(BIN).o
	$(LD) $^ $(LDFLAGS) -o $@

$(LIB)/%.o: $(SRC)/%.cpp
	$(CPP) $(CXXFLAGS) -o $@ -c $<

.PHONY: clean cleanest tidy

clean:
	rm -f $(LIB)/*.o
	rm -f $(BIN)