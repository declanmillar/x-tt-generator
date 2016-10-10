# Fortran makefile
# Declan Millar <d.millar@soton.ac.uk>

UNAME_S := $(shell uname -s)
HOSTNAME := $(shell hostname)

# executable
BIN = zprime

# object files
OBJ = rambo.o vegas.o cteq61pdf.o ct14pdf.o mrs99.o configuration.o lhef.o modelling.o scattering.o alpha_EWNG.o helas.o gg_tt.o qq_tt.o qq_ff.o gg_tt_bbeevv.o qq_tt_bbeevv_qcd.o qq_tt_bbeevv.o gg_bbemuvevm.o qq_bbemuvevm.o uu_bbemuvevm.o dd_bbemuvevm.o rangen.o zprime.o

# gg_bbtatavtvt.o qq_bbtatavtvt.o uu_bbtatavtvt.o dd_bbtatavtvt.o

# Source directory
SRC = Source

# Object directory
LIB = Library

# Output directory
OUT = Binary

ROOTLIBS = -lGui -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -lm -ldl -lTreePlayer

FLAGS = -g -ffree-form -fdefault-real-8 -fdefault-double-8 -std=gnu -J$(LIB) -ffpe-trap=invalid,zero,overflow,underflow,denormal -fmax-errors=0

# Compiler
# ifeq ($(UNAME_S),Linux)
# 	F = ifort
# endif
# ifeq ($(UNAME_S),Darwin)
# 	F = gfortran
# endif

F = gfortran

# Flags
ifeq ($(UNAME_S),Linux)
	FFLAGS = $(FLAGS)
	LFLAGS = -L$(LIB) -lRootTuple -L/afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.08/x86_64-slc6-gcc48-opt/root/lib $(ROOTLIBS)
endif
ifeq ($(HOSTNAME),cyan03)
	FFLAGS = $(FLAGS)
	LFLAGS = -L$(LIB) -lRootTuple -L/local/software/cern/root_v6.06.06/lib $(ROOTLIBS)
endif
ifeq ($(UNAME_S),Darwin)
	# add roottuple libraries os x
	FFLAGS = $(FLAGS)
	LFLAGS = -L$(LIB) -lRootTuple -L/usr/local/Cellar/root6/6.06.08/lib/root $(ROOTLIBS)
endif

# Compile all files ending in .f in SRC
$(LIB)/%.o: $(SRC)/%.f
	$(F) $(FFLAGS) -c -o  $@ $<

# Link mainfile and all processes
$(OUT)/$(BIN): $(patsubst %, $(LIB)/%, $(OBJ))
	$(F) $(FFLAGS) $(LFLAGS) -o $@ $^

# Clean up
clean:
	rm -f $(LIB)/*.o $(LIB)/*.mod $(OUT)/$(BIN)
