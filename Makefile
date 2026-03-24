# define
CXX = g++ -O3 -march=native -ffast-math -funroll-loops -ftree-vectorize
# CXX = g++-14
USE_OPENMP = 1

CXXFLAGS = -Wall -std=c++14
INCLUDES = -Isrc

# if USE_OPENMP=1，add -fopenmp
ifeq ($(USE_OPENMP),1)
    CXXFLAGS += -fopenmp
endif

MAIN_DIR = ./src
OBJ_DIR = ./obj
BIN_DIR = ./bin
TARGET = $(BIN_DIR)/gaukuk.sim

# obtain source files (.cpp files)
SRCS = $(MAIN_DIR)/test.cpp \
	   $(MAIN_DIR)/sim.cpp \
	   $(MAIN_DIR)/eos/adiabatic.cpp \
	   $(MAIN_DIR)/utils/read_config.cpp 

# create object files (.o 文件)
OBJS = $(SRCS:$(MAIN_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# target
all: $(TARGET)

# compile rules
$(OBJ_DIR)/%.o: $(MAIN_DIR)/%.cpp | $(OBJ_DIR)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@


# chain rule
$(TARGET): $(OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

# clean 
clean:
	rm -rf $(OBJ_DIR)/*.o $(TARGET)

# create obj directory (if not exist)
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# create bin directory (if not exist)
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

.PHONY: all clean