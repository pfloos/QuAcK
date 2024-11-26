NVCC = nvcc
NVFLAGS = -O2 --compiler-options '-O2 -Wall -fPIC'

CC = gcc
CFLAGS = -O2 -Wall -g -fPIC

FC = gfortran
FFLAGS = -O2 -Wall -g -fPIC

SRC_DIR = src
INC_DIR = include
BLD_DIR = build
$(shell mkdir -p $(BLD_DIR))

CU_SRC = $(wildcard $(SRC_DIR)/*.cu)
CU_OBJ = $(CU_SRC:$(SRC_DIR)/%.cu=$(BLD_DIR)/%.o)

C_SRC = $(wildcard $(SRC_DIR)/*.c)
C_OBJ = $(C_SRC:$(SRC_DIR)/%.c=$(BLD_DIR)/%.o)

F_SRC = #$(SRC_DIR)/cu_quack_module.f90
F_OBJ = #$(BLD_DIR)/cu_quack_module.o

OUTPUT_LIB = $(BLD_DIR)/libcuquack.so

all: $(OUTPUT_LIB)

$(OUTPUT_LIB): $(CU_OBJ) $(C_OBJ) $(F_OBJ)
	$(CC) -shared -o $(OUTPUT_LIB) $(CU_OBJ) $(C_OBJ) $(F_OBJ)

$(BLD_DIR)/%.o: $(SRC_DIR)/%.cu
	$(NVCC) $(NVFLAGS) -c -o $@ $< -I$(INC_DIR)

$(BLD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $< -I$(INC_DIR)

$(F_OBJ): $(F_SRC)
	$(FC) $(FFLAGS) -c -o $@ $< -J$(BLD_DIR)

.PHONY: clean
clean:
	rm -f $(BLD_DIR)/*


