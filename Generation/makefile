# Fortran makefile
# Declan Millar <d.millar@soton.ac.uk>

# executable
BIN = zprime

# object files
OBJ = vegas.o cteq61pdf.o mathematics.o configuration.o modelling.o kinematics.o class_histogram.o class_histogram2d.o distributions.o  alpha_EWNG.o helas.o dsigma.o ggbbffff_qcd.o ggff_qcd.o mrs99.o qqbbffff_ewp.o qqbbffff_qcd.o qqff_ewp.o qqff_qcd.o rangen.o zprime.o

# Source directory
SRC = Source

# Object directory
LIB = Library

# Output directory
OUT = Binary

# Compiler
F = gfortran

# Flags
FFLAGS = -g -ffree-form -fdefault-real-8 -fdefault-double-8 -std=gnu -J$(LIB) -ffpe-trap=invalid,zero,overflow,underflow -fmax-errors=0
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	# add roottuple libraries qmulpc007
	LFLAGS = -L/usr/local/lib -lRootTuple -L/afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.28/x86_64-slc6-gcc48-opt/root/lib -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lTreePlayer
endif
ifeq ($(UNAME_S),Darwin)
	# add roottuple libraries os x
	LFLAGS = -L/usr/lib -lRootTuple -L/usr/local/Cellar/root/5.34.26/lib/root -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -lm -ldl -lTreePlayer -Wl,-stack_size,0x40000000
endif

# Compile all files ending in .f in SRC
$(LIB)/%.o: $(SRC)/%.f 
	$(F) $(FFLAGS) -c -o  $@ $<

# Link mainfile and all processes
$(OUT)/$(BIN): $(patsubst %, $(LIB)/%, $(OBJ))
	$(F) $(FFLAGS) $(LFLAGS) -o $@ $^

# Clean up
clean:
	rm -f $(LIB)/*.o $(OUT)/$(BIN) 
