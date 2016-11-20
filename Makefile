# fortran makefile for zprime-top-generator
# needs RootTuple libraries
# author: Declan Millar <declan.millar@cern.ch>

F = gfortran
SRC = src
LIB = lib
OUT = bin
BIN = generator

OBJ = vamp_kinds.o progress.o exceptions.o vamp_stat.o utils.o divisions.o histograms.o iso_fortran_env_stub.o linalg.o products.o specfun.o tao52_random_numbers.o tao_random_numbers.o vamp.o rambo.o cteq61pdf.o ct14pdf.o mrs99.o configuration.o lhef.o modelling.o scattering.o alpha_EWNG.o helas.o gg_tt.o qq_tt.o qq_ff.o gg_tt_bbeevv.o qq_tt_bbeevv_qcd.o qq_tt_bbeevv.o gg_bbemuvevm.o qq_bbemuvevm.o uu_bbemuvevm.o dd_bbemuvevm.o rangen.o generator.o


ROOTLIBS = -lGui -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -lm -ldl -lTreePlayer
CFLAGS = -J$(LIB) -ffree-form -ffpe-trap=invalid,zero,overflow,underflow,denormal 

OS := $(shell uname)
HOST := $(shell hostname)
# ifeq ($(HOST),Linux)

ifneq (,$(findstring lxplus,$(HOST)))
	LFLAGS = -L$(LIB) -lRootTuple -L/afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.08/x86_64-slc6-gcc49-opt/root/lib $(ROOTLIBS)
else ifneq (,$(findstring cyan,$(HOST)))
	LFLAGS = -L$(LIB) -lRootTuple -L/local/software/cern/root_v6.06.06/lib $(ROOTLIBS)
else ifneq (,$(findstring heppc,$(HOST)))
	LFLAGS = -L$(LIB) -lRootTuple -L/opt/root/lib $(ROOTLIBS)
else ifneq (,$(findstring Sunder,$(HOST)))
	LFLAGS = -L$(LIB) -lRootTuple -L/usr/local/Cellar/root6/6.06.08/lib/root $(ROOTLIBS)
else
test:
	$(info *** I do not know where to find the ROOT libraries! ***)
endif

# Compile all files ending in .f in SRC
$(LIB)/%.o: $(SRC)/%.f
	$(F) $(CFLAGS) -c -o  $@ $<

# Link main file and all processes
$(OUT)/$(BIN): $(patsubst %, $(LIB)/%, $(OBJ))
	$(F) $(LFLAGS) -o $@ $^

clean:
	rm -f $(LIB)/*.o $(LIB)/*.mod $(OUT)/$(BIN)
