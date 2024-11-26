NVCC = nvcc
NFLAGS = -O2 --compiler-options '-fPIC'
NDFLAGS = --shared

CC = gcc
CFLAGS = -fPIC -O2 -Wall -g

FC = gfortran
FFLAGS = -O2 -Wall -g

SRC_DIR = src
INC_DIR = include

BIN_DIR = bin
BLD_DIR = build
$(shell mkdir -p $(BIN_DIR))
$(shell mkdir -p $(BLD_DIR))

CU_SRC = $(wildcard $(SRC_DIR)/*.cu)
CU_OBJ = $(CU_SRC:$(SRC_DIR)/%.cu=$(BLD_DIR)/%.o)

C_SRC = $(SRC_DIR)/ph_drpa.c
C_OBJ = $(BLD_DIR)/ph_drpa.o

F_SRC = $(SRC_DIR)/cu_quack_module.f90
F_OBJ = $(BLD_DIR)/cu_quack_module.f90

MAIN_SRC = $(SRC_DIR)/cu_quack.f90
MAIN_OBJ = $(BLD_DIR)/cu_quack.o

OUTPUT_LIB = $(BLD_DIR)/libcuquack.so

CUDA_LIBS = -lcudart -lcublas


all: $(OUTPUT_LIB)

$(OUTPUT_LIB): $(CU_OBJ) $(C_OBJ)
	$(NVCC) $(NFLAGS) $(NLDFLAGS) $^ -o $@ $(CUDA_LIBS) -I$(INC_DIR)

$(BLD_DIR)/%.o: $(SRC_DIR)/%.cu
	$(NVCC) $(NFLAGS) -c $< -o $@ -I$(INC_DIR)

$(C_OBJ): $(C_SRC)
	@for src in $(C_SRC); do \
		obj=$(BLD_DIR)/$$(basename $${src} .c).o; \
		echo "$(CC) $(CFLAGS) -c $$src -o $$obj -I$(INC_DIR)"; \
		$(CC) $(CFLAGS) -c $$src -o $$obj -I$(INC_DIR); \
	done

$(F_OBJ): $(F_SRC)
	$(FC) $(FFLAGS) -c $< -o $@ -J$(BLD_DIR)

$(MAIN_OBJ): $(MAIN_SRC)
	$(FC) $(FFLAGS) -c $< -o $@ -J$(BLD_DIR)



.PHONY: clean
clean:
	rm -f $(BLD_DIR)/*.o $(BLD_DIR)/*.so $(BLD_DIR)/*.mod $(BIN_DIR)/*


