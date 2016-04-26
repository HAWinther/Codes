#####################################################
# Hans A. Winther (2015) (hans.a.winther@gmail.com)
#####################################################

SHELL := /bin/bash

# Set compiler
CC = g++ -std=c++11

C = -O3 $(OPTIONS)

TARGETS := qubicspline
all: $(TARGETS)

# OBJECT FILES
OBJS = Spline.o main.o

# DEPENDENCIES
main.o		: main.cpp Spline.cpp

# HEADERS
HEADERS = Spline.h

qubicspline: $(OBJS)
	${CC} -o $@ $^ $C $I $L

%.o: %.cpp $(HEADERS)
	${CC}  -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o

