SHELL := /bin/bash

C = -O3 $(OPTIONS)

TARGETS := readwriteramses_test
all: $(TARGETS)

# OBJECT FILES
OBJS = main.o

# DEPENDENCIES
main.o		: main.cpp read_write_ramses.h

# HEADERS
HEADERS = read_write_ramses.h

readwriteramses_test: $(OBJS)
	g++-mp-4.8 -o $@ $^ $C $I $L

%.o: %.cpp $(HEADERS)
	g++-mp-4.8  -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o

