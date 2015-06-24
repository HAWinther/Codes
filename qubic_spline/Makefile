SHELL := /bin/bash

# Define the type
type=DOUBLE

C = -O3 $(OPTIONS) -D$(type)

TARGETS := qubicspline
all: $(TARGETS)

# OBJECT FILES
OBJS = main.o

# DEPENDENCIES
main.o		: main.cpp spline.h

# HEADERS
HEADERS = spline.h

qubicspline: $(OBJS)
	g++-mp-4.8 -o $@ $^ $C $I $L

%.o: %.cpp $(HEADERS)
	g++-mp-4.8  -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o

