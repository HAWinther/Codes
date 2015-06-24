SHELL := /bin/bash

# Define the type
type=DOUBLE

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
	g++-mp-4.8 -o $@ $^ $C $I $L

%.o: %.cpp $(HEADERS)
	g++-mp-4.8  -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o

