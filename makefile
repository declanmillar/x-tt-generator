# ==============================================================================
# Fortran makefile
#
# Authors: Declan Millar
# ------------------------------------------------------------------------------
# Definitions

# executable
BIN = zprime

# object files
OBJ = vegas.o cteq61pdf.o configuration.o kinematics.o distributions.o  alpha_EWNG.o coupZp.o helas.o differential_cross_section.o ggbbffff_qcd.o ggff_qcd.o initialise_madgraph.o mrs99.o qqbbffff_ewp.o qqbbffff_qcd.o qqff_ewp.o qqff_qcd.o rangen.o width_zprime.o zprime.o

# Source directory
SRC = Source

# Object directory
LIB = Library

# Output directory
OUT = Binary

# Compiler
F = gfortran

# Flags
FFLAGS = -g -ffree-form -fdefault-real-8 -fdefault-double-8 -std=gnu -ffpe-trap=invalid,zero,overflow,underflow -fmax-errors=0 -J$(LIB)
# -Wall
# ------------------------------------------------------------------------------
# Commands

# Compile all files ending in .f in SRC
$(LIB)/%.o: $(SRC)/%.f 
	$(F) $(FFLAGS) -c -o  $@ $<

# Link mainfile and all processes
$(OUT)/$(BIN): $(patsubst %,$(LIB)/%, $(OBJ))
	$(F) $(FFLAGS) -o $@ $^

# Clean up
clean:
	rm -f $(LIB)/*.o $(OUT)/$(BIN)

# ============================================================================== 
