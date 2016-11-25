# fortran makefile for zprime-top-generator
# requires: RootTuple
# author: Declan Millar <declan.millar@cern.ch>

F = gfortran
SRC = src
LIB = lib
OUT = bin
BIN = generator

OBJ = vamp_kinds.o progress.o exceptions.o vamp_stat.o utils.o divisions.o histograms.o iso_fortran_env_stub.o linalg.o products.o specfun.o tao52_random_numbers.o tao_random_numbers.o vamp.o rambo.o cteq61pdf.o ct14pdf.o mrs99.o configuration.o lhef.o helas.o modelling.o scattering.o alpha_EWNG.o gg_tt.o qq_tt.o qq_ff.o gg_tt_bbeevv.o qq_tt_bbeevv_qcd.o qq_tt_bbeevv.o gg_bbemuvevm.o qq_bbemuvevm.o uu_bbemuvevm.o dd_bbemuvevm.o generator.o

FFLAGS = -J$(LIB) -std=f2008 -ffpe-trap=invalid,zero,overflow,underflow,denormal 
LFLAGS = -L$(LIB) -lRootTuple $(shell root-config --libs)

$(LIB)/%.o: $(SRC)/%.f03
	$(F) $(FFLAGS) -c -o  $@ $<

$(OUT)/$(BIN): $(patsubst %, $(LIB)/%, $(OBJ))
	$(F) $(LFLAGS) -o $@ $^

clean:
	rm -f $(LIB)/*.o $(LIB)/*.mod $(OUT)/$(BIN)
