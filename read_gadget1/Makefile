SHELL := /bin/bash

# Set compiler
# CC = g++-mp-4.8

C = -O3

TARGETS := readgadgetfile
all: $(TARGETS)

# OBJECT FILES
OBJS = main.o

# DEPENDENCIES
main.o		: main.cpp read_gadget.h

# HEADERS
HEADERS = read_gadget.h

readgadgetfile: $(OBJS)
	${CC} -o $@ $^ $C $I $L

%.o: %.cpp $(HEADERS)
	${CC}  -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o

