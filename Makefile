# fortran makefile for zprime-top-generator
# requires: RootTuple
# author: Declan Millar <declan.millar@cern.ch>

# HOSTNAME := $(shell hostname)
# ifneq (,$(findstring lxplus, $(HOSTNAME)))
    # F = ifort
# else
F = gfortran
F = mpifort
# endif

SRC = src
LIB = lib
OUT = bin
BIN = generator
TUPLELIB = ../root-tuple/lib

OBJ = vamp_kinds.o runtime.o mpi90.o constants.o specfun.o coordinates.o progress.o exceptions.o vamp_stat.o utils.o divisions.o histograms.o iso_fortran_env_stub.o linalg.o products.o specfun.o tao52_random_numbers.o tao_random_numbers.o vamp.o vampi.o vamp_tests.o vamp_test0.o rambo.o cteq61pdf.o ct14pdf.o mrs99.o configuration.o lhef.o helas.o modelling.o alpha_EWNG.o gg_tt.o qq_tt.o qq_ff.o tt_bbeevv.o gg_bbemuvevm.o qq_bbemuvevm.o uu_bbemuvevm.o dd_bbemuvevm.o scattering.o generator.o

# ifneq (,$(findstring lxplus, $(HOSTNAME)))
    # FFLAGS = -warn all
# else
FFLAGS = -J$(LIB) -std=f2008 -ffpe-trap=invalid,zero,overflow,underflow,denormal
# endif

LFLAGS = -L$(LIB) -L$(TUPLELIB) -lRootTuple $(shell root-config --libs)

$(LIB)/%.o: $(SRC)/%.f90
	$(F) $(FFLAGS) -c -o $@ $<

$(OUT)/$(BIN): $(patsubst %, $(LIB)/%, $(OBJ))
	$(F) $(LFLAGS) -o $@ $^

clean:
	rm -f $(LIB)/*.o $(LIB)/*.mod $(OUT)/$(BIN)
