# define
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
    # macOS → use clang
    CXX = clang++
    OPENMP_FLAG = -Xpreprocessor -fopenmp
    OPENMP_LIB  = -lomp
    OPENMP_INC  = -I/opt/homebrew/opt/libomp/include
    OPENMP_LIBPATH = -L/opt/homebrew/opt/libomp/lib
else
    # Linux → use g++
    CXX = g++
    OPENMP_FLAG = -fopenmp
    OPENMP_LIB  =
    OPENMP_INC  =
    OPENMP_LIBPATH =
endif 
USE_OPENMP = 1

CXXFLAGS = -O3 -march=native -ffast-math -funroll-loops -std=c++14 -Wall \
#            -Rpass=loop-vectorize
INCLUDES = -Isrc
LDFLAGS =

# if USE_OPENMP=1，add -fopenmp
ifeq ($(USE_OPENMP),1)
    CXXFLAGS += $(OPENMP_FLAG) $(OPENMP_INC)
    LDFLAGS  += $(OPENMP_LIBPATH) $(OPENMP_LIB)
endif

# paths 
MAIN_DIR = ./src
OBJ_DIR = ./obj
BIN_DIR = ./bin
TARGET = $(BIN_DIR)/gaukuk.sim

# obtain source files (.cpp files)
SRCS = $(MAIN_DIR)/test.cpp \
	   $(MAIN_DIR)/sim.cpp \
	   $(MAIN_DIR)/eos/adiabatic.cpp \
	   $(MAIN_DIR)/utils/read_config.cpp 
# 	   $(MAIN_DIR)/grid/

# create object files (.o 文件)
OBJS = $(SRCS:$(MAIN_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# target
all: $(TARGET)

# compile rules
$(OBJ_DIR)/%.o: $(MAIN_DIR)/%.cpp | $(OBJ_DIR)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# chain rule
$(TARGET): $(OBJS) | $(BIN_DIR)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(OBJS) $(LDFLAGS) $(OPENMP_FLAG) -o $(TARGET)

# clean 
clean:
	rm -rf $(OBJ_DIR)/* $(TARGET)

# create obj directory (if not exist)
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# create bin directory (if not exist)
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

.PHONY: all clean