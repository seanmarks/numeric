
CXX := g++

# Basic flags
FLAGS := -g -std=c++11

# GCC-specific OpenMP flags
#GCC_FLAGS := -fopenmp -fopenmp-simd
GCC_FLAGS := -fopenmp-simd

# -fopt-info-vec-missed
#    Note: missed/optimized/all
#    Older: -ftree-vectorizer-verbose=1
#GCC_FLAGS += -ftree-vectorize 
#GCC_FLAGS += -fopt-info-vec-all
#GCC_FLAGS += -fopt-info-vec-optimized
#GCC_FLAGS += -fopt-info-vec-missed

# Optimizations
#OPT_FLAGS := -O2
OPT_FLAGS := -O3
OPT_FLAGS += -march=native

# Combine flags
CXXFLAGS := $(FLAGS) $(GCC_FLAGS) $(OPT_FLAGS) 

PROJECT := test

.PHONY: all $(PROJECT) clean
all: $(PROJECT)

$(PROJECT):
	$(CXX) -o $@ $(CXXFLAGS) main.cpp System.cpp

clean:
	rm $(PROJECT)

#$(CXX) -S -o main.s $(CXXFLAGS) main.cpp
#$(CXX) -S -o main.s -fverbose-asm $(CXXFLAGS) main.cpp
