# Makefile

PETSC_ARCH=/chalmers/sw/sup64/petsc-3.7.4

include ${PETSC_ARCH}/lib/petsc/conf/variables
include ${PETSC_ARCH}/lib/petsc/conf/rules

CXX=g++
CXXFLAGS=-Wall -Wextra -g -O0 -c -Iinclude -I${PETSC_ARCH}/include
LD=g++
LFLAGS=

OBJECTS=Main.o CG.o Create.o DiscretePoisson2D.o GaussSeidel.o Jacobi.o  PCG.o Solver.o SOR.o
Run=Main

all: $(Run)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

$(Run): $(OBJECTS)
	$(LD) $(LFLAGS) $(OBJECTS) $(PETSC_LIB) -o $@
