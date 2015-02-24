# ==============================================================================
# Fortran makefile
#
# Authors: Stefano Moretti, Declan Millar
# ------------------------------------------------------------------------------
# Definitions

# executable
BIN = zprime

# Source directory
SRC = .

# Output directory
OUT = ..

# Object directory
OBJ = .

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
zprime: *.o makefile
	$(F) -o $(OUT)/$(BIN) $(FFLAGS) *.o

# Clean up
clean:
	rm -f $(SRC)/*.o

# ============================================================================== 
