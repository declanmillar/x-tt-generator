# ==============================================================================
# Fortran makefile
#
# Authors: Declan Millar
# ------------------------------------------------------------------------------
# Definitions

# executable
BIN = zprime

# object files
OBJ = Cteq61Pdf.o alpha_EWNG.o coupZp.o dhelas_all.o dxsec.o ggbbffff_qcd.o ggff_qcd.o initialise_madgraph.o mrs99.o qqbbffff_ewp.o qqbbffff_qcd.o qqff_ewp.o qqff_qcd.o rangen.o ve_dist.o widthZp.o zprime.o

# Source directory
SRC = Source

# Object directory
LIB = Library

# Output directory
OUT = Binary

# Compiler
F = gfortran

# Flags
FFLAGS = -g -ffpe-trap=invalid,zero,overflow,underflow,denormal
# -Wall
# ------------------------------------------------------------------------------
# Commands

# Compile all files ending in .f in SRC
$(LIB)/%.o: $(SRC)/%.f 
	$(F) $(FFLAGS) -c -o $@ $<

# Link mainfile and all processes
$(OUT)/$(BIN): $(patsubst %,$(LIB)/%, $(OBJ))
	$(F) $(FFLAGS) -o $@ $^

# Clean up
clean:
	rm -f $(LIB)/*.o $(OUT)/$(BIN)

# ============================================================================== 
