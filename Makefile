# Compiler
FC = ifx

# Directories
SRC_DIR = src
OBJ_DIR = bin\obj
MOD_DIR = bin\mods
BIN_DIR = bin

# Target executable
TARGET = lbm.exe

# Compiler flags
FFLAGS = /nologo /module:$(MOD_DIR)

# Sources (dependency order matters for Fortran modules)
SRC = \
	$(SRC_DIR)/lbm_consts.f90 \
	$(SRC_DIR)/lbm_solver.f90 \
	$(SRC_DIR)/lbm_bc.f90 \
	$(SRC_DIR)/lbm_utils.f90 \
	$(SRC_DIR)/lbm_simulation.f90 \
	$(SRC_DIR)/lbm_config.f90 \
	$(SRC_DIR)/lbm_swd.f90 \
	$(SRC_DIR)/lbm_dwd.f90 \
	$(SRC_DIR)/lbm_io.f90 \
	$(SRC_DIR)/lbm.f90

OBJ = $(SRC:$(SRC_DIR)/%.f90=$(OBJ_DIR)/%.obj)

# Default target
all: dirs $(TARGET)

# Create directories
dirs:
	@if not exist $(OBJ_DIR) mkdir $(OBJ_DIR)
	@if not exist $(MOD_DIR) mkdir $(MOD_DIR)

# Linking
$(TARGET): $(OBJ)
	$(FC) /nologo $(OBJ) /exe:$@

# Compilation
$(OBJ_DIR)/%.obj: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) /c $< /object:$@

# Remove build files
clean:
	del /Q $(OBJ_DIR)\*.obj
	del /Q $(MOD_DIR)\*.mod
	del /Q $(TARGET)

run: all
	$(TARGET)