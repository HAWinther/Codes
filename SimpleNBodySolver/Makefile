SHELL := /bin/bash

##################################################

# Options
OPTIONS = -O3 -std=c++11

# Path to FFTW3 includes
I  = -I$(HOME)/local/include
L  = -L$(HOME)/local/lib -lfftw3  # -lfftw3_threads 

# LUA includes
I += -I/Users/winther/local/lua-5.2.3/install/include
L += -L/Users/winther/local/lua-5.2.3/install/lib  -llua -ldl

# OpenMP?
OPTIONS += -DOMP -fopenmp        # -DOMP_FFT

# TSC density assignment
OPTIONS += -DTSC

# f(R) gravity [if this is not defined the MG solver solves for gravity potential]
OPTIONS += -DSYMMETRON

# Set C++ compiler
CC = g++-mp-4.8
C = $(OPTIONS)

##################################################

TARGETS := nbodysolver
all: $(TARGETS)
OBJS = Spline.o main.o
main.o		: main.cpp
HEADERS = global.h Spline.h multigrid_solver.h

nbodysolver: $(OBJS)
	${CC} -o $@ $^ $C $I $L

%.o: %.cpp $(HEADERS)
	${CC}  -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o

