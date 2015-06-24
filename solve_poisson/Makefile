SHELL := /bin/bash

# Set the type
type=DOUBLE

# Includes for FFTW3
L = -L/Users/winther/local/lib
I = -I/Users/winther/local/include -lm -lfftw3

# For my MAC, choose correct compiler:
# CC = g++-mp-4.8 

C = -O3 $(OPTIONS) -D$(type)

TARGETS := poissonsolver
all: $(TARGETS)

# OBJECT FILES
OBJS = main.o

# DEPENDENCIES
main.o		: solve_poisson.h

# HEADERS
HEADERS = solve_poisson.h

poissonsolver: $(OBJS)
	$(CC) -o $@ $^ $C $I $L $(SFML)

%.o: %.cpp $(HEADERS)
	$(CC) -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o

