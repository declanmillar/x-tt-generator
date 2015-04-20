# Fortran makefile
# Declan Millar <d.millar@soton.ac.uk>

# executable
BIN = zprime

# object files
OBJ = vegas.o cteq61pdf.o mathematics.o configuration.o modelling.o kinematics.o class_histogram.o class_histogram2d.o distributions.o  alpha_EWNG.o helas.o differential_cross_section.o ggbbffff_qcd.o ggff_qcd.o mrs99.o qqbbffff_ewp.o qqbbffff_qcd.o qqff_ewp.o qqff_qcd.o rangen.o zprime.o

# Source directory
SRC = Source

# Object directory
LIB = Library

# Output directory
OUT = Binary

# Compiler
F = gfortran

# Flags
FFLAGS = -g -ffree-form -fdefault-real-8 -fdefault-double-8 -std=gnu -J$(LIB) -ffpe-trap=invalid,zero,overflow,underflow -fmax-errors=0 
LFLAGS = -Wl,-stack_size,0x40000000

# Compile all files ending in .f in SRC
$(LIB)/%.o: $(SRC)/%.f 
	$(F) $(FFLAGS) -c -o  $@ $<

# Link mainfile and all processes
$(OUT)/$(BIN): $(patsubst %, $(LIB)/%, $(OBJ))
	$(F) $(FFLAGS) $(LFLAGS) -o $@ $^

# Clean up
clean:
	rm -f $(LIB)/*.o $(OUT)/$(BIN) 
