#####################################################
# Hans A. Winther (2015) (hans.a.winther@gmail.com)
#####################################################

SHELL := /bin/bash

# Set compiler
CC = g++-mp-4.8 -std=c++11

C = -O3

TARGETS := mglinearpert
all: $(TARGETS)

# OBJECT FILES
OBJS = Spline.o main.o

# DEPENDENCIES
Spline.o  : Spline.cpp
main.o		: main.cpp

# HEADERS
HEADERS = PofkModifiedGravity.h Spline.h OdeSolver.h

mglinearpert: $(OBJS)
	${CC} -o $@ $^ $C $I $L

%.o: %.cpp $(HEADERS)
	${CC} -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o

