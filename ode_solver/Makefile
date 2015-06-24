SHELL := /bin/bash

# Define the type
type=DOUBLE

# Set compiler
# CC = g++-mp-4.8

C = -O3 $(OPTIONS) -D$(type)

TARGETS := odesolver
all: $(TARGETS)

# OBJECT FILES
OBJS = main.o

# DEPENDENCIES
main.o		: main.cpp ode_solver.h

# HEADERS
HEADERS = ode_solver.h

odesolver: $(OBJS)
	${CC} -o $@ $^ $C $I $L

%.o: %.cpp $(HEADERS)
	${CC}  -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o

