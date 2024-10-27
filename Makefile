# Compiler and flags
FC = gfortran
CFLAGS = -O3 -fopenmp -J$(BIN_DIR)  # -J specifies the directory for .mod files
LDFLAGS = -fopenmp

# Directories
SRC = main.f90   # Replace or add additional source files
BIN_DIR = bin
OBJ = $(patsubst %.f90,$(BIN_DIR)/%.o,$(SRC))  # Place .o files in bin
TARGET = $(BIN_DIR)/Ternary_Contact_Angle_On_Flat_Surface  # Target executable in bin directory

# Ensure bin directory exists and build the target with ulimit set
$(TARGET): $(OBJ) | $(BIN_DIR)
	ulimit -s unlimited && $(FC) $(LDFLAGS) -o $(TARGET) $(OBJ)

# Rule to create bin directory if it doesn't exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Compile each .f90 file to .o in bin
$(BIN_DIR)/%.o: %.f90 | $(BIN_DIR)
	$(FC) $(CFLAGS) -c $< -o $@
	
# Run target with specified number of threads
run: $(TARGET)
	OMP_NUM_THREADS=$(NUM_THREADS) ./$(TARGET)

# Clean target
.PHONY: clean
clean:
	rm -rf $(BIN_DIR) *.csv *.vtk *.vtr
	
	# Default target
.DEFAULT_GOAL := run
