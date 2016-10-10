# Fortran makefile for zprime-top-generator
# Author: Declan Millar <declan.millar@cern.ch>

F = gfortran
SRC = Source
LIB = Library
OUT = Binary
BIN = zprime
OBJ = rambo.o vegas.o cteq61pdf.o ct14pdf.o mrs99.o configuration.o lhef.o modelling.o scattering.o alpha_EWNG.o helas.o gg_tt.o qq_tt.o qq_ff.o gg_tt_bbeevv.o qq_tt_bbeevv_qcd.o qq_tt_bbeevv.o gg_bbemuvevm.o qq_bbemuvevm.o uu_bbemuvevm.o dd_bbemuvevm.o rangen.o zprime.o

# gg_bbtatavtvt.o qq_bbtatavtvt.o uu_bbtatavtvt.o dd_bbtatavtvt.o

ROOTLIBS = -lGui -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -lm -ldl -lTreePlayer
CFLAGS = -g -ffree-form -fdefault-real-8 -fdefault-double-8 -std=gnu -J$(LIB) -ffpe-trap=invalid,zero,overflow,underflow,denormal -fmax-errors=0

OS := $(shell uname)
HOST := $(shell hostname)
ifeq ($(OS),Linux)
	LFLAGS = -L$(LIB) -lRootTuple -L/afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.08/x86_64-slc6-gcc48-opt/root/lib $(ROOTLIBS)
endif
ifeq ($(HOST),cyan03)
	LFLAGS = -L$(LIB) -lRootTuple -L/local/software/cern/root_v6.06.06/lib $(ROOTLIBS)
endif
ifeq ($(OS),Darwin)
	LFLAGS = -L$(LIB) -lRootTuple -L/usr/local/Cellar/root6/6.06.08/lib/root $(ROOTLIBS)
endif

# Compile all files ending in .f in SRC
$(LIB)/%.o: $(SRC)/%.f
	$(F) $(CFLAGS) -c -o  $@ $<

# Link main file and all processes
$(OUT)/$(BIN): $(patsubst %, $(LIB)/%, $(OBJ))
	$(F) $(LFLAGS) -o $@ $^

clean:
	rm -f $(LIB)/*.o $(LIB)/*.mod $(OUT)/$(BIN)
