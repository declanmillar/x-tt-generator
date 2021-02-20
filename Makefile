# Author: Declan Millar

OBJ = vamp_kinds.o runtime.o constants.o specfun.o coordinates.o progress.o exceptions.o vamp_stat.o utils.o iso_fortran_env_stub.o divisions.o histograms.o linalg.o products.o specfun.o tao52_random_numbers.o tao_random_numbers.o vamp.o vamp_tests.o vamp_test0.o rambo.o cteq61pdf.o ct14pdf.o mrs99.o configuration.o lhef.o helas.o modelling.o alpha_EWNG.o gg_tt.o qq_tt.o qq_ff.o tt_bbeevv.o gg_bbemuvevm.o qq_bbemuvevm.o uu_bbemuvevm.o dd_bbemuvevm.o scattering.o generator.o
BIN = generator
SRC_DIR = src/
LIB_DIR = lib/
BIN_DIR = bin/

IFORT := $(shell command -v ifort 2> /dev/null)
GFORTRAN := $(shell command -v gfortran 2> /dev/null)

ifdef IFORT
	F = ifort
	FFLAGS = -stand f08 -parallel -qopenmp -O3 -module $(LIB_DIR)
	LFLAGS = -parallel -qopenmp -L$(LIB_DIR)
else ifdef GFORTRAN
	F = gfortran
	FFLAGS = -std=f2008 -ffpe-trap=invalid,zero,overflow,underflow,denormal -pedantic -Wall -J$(LIB_DIR)
	LFLAGS = -L$(LIB_DIR)
else
	$(error "Found neither ifort nor gfortran. Install ifort or gfortran.")
endif

$(LIB_DIR)%.o: $(SRC_DIR)%.f90
	@mkdir -p $(@D)
	$(F) $(FFLAGS) -c -o $@ $<

$(BIN_DIR)$(BIN): $(patsubst %, $(LIB_DIR)%, $(OBJ))
	@mkdir -p $(@D)	
	$(F) $(LFLAGS) -o $@ $^
	
.PHONY: clean
clean:
	rm -rf $(LIB_DIR) $(BIN_DIR)
