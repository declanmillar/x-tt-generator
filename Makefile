# fortran makefile for perigee
# requires: root-tuple as a sibling with perigee
# author: Declan Millar <declan.millar@cern.ch>

BIN = generator
SRC = src
LIB = lib
OUT = bin
TUPLELIB = ../root-tuple/lib

HOSTNAME := $(shell hostname)
ifeq ($(HOSTNAME), Lorkhan)
	F = gfortran
	L = gfortran
	FFLAGS = -J$(LIB) -std=f2008 -ffpe-trap=invalid,zero,overflow,underflow,denormal
else
	F = ifort
	L = ifort
	FFLAGS = -module $(LIB) -parallel -openmp
endif

OBJ = vamp_kinds.o runtime.o constants.o specfun.o coordinates.o progress.o exceptions.o vamp_stat.o utils.o iso_fortran_env_stub.o divisions.o histograms.o linalg.o products.o specfun.o tao52_random_numbers.o tao_random_numbers.o vamp.o vamp_tests.o vamp_test0.o rambo.o cteq61pdf.o ct14pdf.o mrs99.o configuration.o lhef.o helas.o modelling.o alpha_EWNG.o gg_tt.o qq_tt.o qq_ff.o tt_bbeevv.o gg_bbemuvevm.o qq_bbemuvevm.o uu_bbemuvevm.o dd_bbemuvevm.o scattering.o generator.o

LFLAGS = -L$(LIB) -L$(TUPLELIB) -lRootTuple

$(LIB)/%.o: $(SRC)/%.f90
	$(F) $(FFLAGS) -c -o $@ $<

$(OUT)/$(BIN): $(patsubst %, $(LIB)/%, $(OBJ))
	$(L) $(LFLAGS) -o $@ $^

clean:
	rm -f $(LIB)/*.o $(LIB)/*.mod $(OUT)/$(BIN)
