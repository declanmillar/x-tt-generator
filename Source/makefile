# ==============================================================================
# Fortran makefile
#
# Authors: Stefano Moretti, Declan Millar
# ------------------------------------------------------------------------------
# Definitions

# executable
BIN = ttbar_BSM

# Source directory
SRC = .

# Output directory
OUT = ..

# Compiler
F = gfortran

# Flags
FFLAGS = -g -ffpe-trap=invalid,zero,overflow,underflow,denormal
# -Wall
# ------------------------------------------------------------------------------
# Commands

# Compile all files ending in .f
.f.o: ; $(F) $(FFLAGS) -c $(SRC)/*.f

# Link mainfile and all processes
ppbbllnn_SM: *.o makefile
	$(F) -o $(OUT)/$(BIN) $(FFLAGS) *.o

# Clean up
clean:
	rm -f $(SRC)/*.o

# ============================================================================== 
