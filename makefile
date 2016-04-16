# Fortran makefile
# Declan Millar <d.millar@soton.ac.uk>

UNAME_S := $(shell uname -s)

# executable
BIN = zprime

# object files
OBJ = rambo.o vegas.o cteq61pdf.o mrs99.o configuration.o modelling.o scattering.o alpha_EWNG.o helas.o gg_tt.o qq_tt.o qq_ff.o gg_tt_bbeevv.o qq_tt_bbeevv_qcd.o qq_tt_bbeevv.o gg_bbemuvevm.o qq_bbemuvevm.o uu_bbemuvevm.o dd_bbemuvevm.o rangen.o zprime.o

# Source directory
SRC = Source

# Object directory
LIB = Library

# Output directory
OUT = Binary

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
	# add roottuple libraries qmulpc007
	# FFLAGS = -g -real_size 64 -double_size 64 -free -module $(LIB)
	FFLAGS = -g -ffree-form -fdefault-real-8 -fdefault-double-8 -std=gnu -J$(LIB) -ffpe-trap=invalid,zero,overflow,underflow,denormal -fmax-errors=0
	LFLAGS = -L/afs/cern.ch/user/d/demillar/.RootTuple -lRootTuple -L/afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.28/x86_64-slc6-gcc48-opt/root/lib -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lTreePlayer
endif
ifeq ($(UNAME_S),Darwin)
	# add roottuple libraries os x
	FFLAGS = -g -ffree-form -fdefault-real-8 -fdefault-double-8 -std=gnu -J$(LIB) -ffpe-trap=invalid,zero,overflow,underflow,denormal -fmax-errors=0
	LFLAGS = -L/usr/local/lib -lRootTuple -L/usr/local/Cellar/root/5.34.34_1/lib/root -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -lm -ldl -lTreePlayer
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
