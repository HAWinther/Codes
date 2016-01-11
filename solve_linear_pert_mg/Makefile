#####################################################
# Hans A. Winther (2015) (hans.a.winther@gmail.com)
#####################################################

SHELL := /bin/bash

# Define the type
type=DOUBLE

# Set compiler
CC = g++-mp-4.8

C = -O3 -D$(type)

TARGETS := mglinearpert
all: $(TARGETS)

# OBJECT FILES
OBJS = main.o

# DEPENDENCIES
main.o		: main.cpp linear_pert_solver.h spline.h ode_solver.h

# HEADERS
HEADERS = linear_pert_solver.h spline.h ode_solver.h

mglinearpert: $(OBJS)
	${CC} -o $@ $^ $C $I $L

%.o: %.cpp $(HEADERS)
	${CC} -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o

