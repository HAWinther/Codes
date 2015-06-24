SHELL := /bin/bash

# Set compiler
CC = g++-mp-4.8

C = -O3

# Set paths
L = -L/Users/winther/local/lua-5.2.3/install/lib  -llua -ldl
I = -I/Users/winther/local/lua-5.2.3/install/include

TARGETS := readluafile
all: $(TARGETS)

# OBJECT FILES
OBJS = main.o

# DEPENDENCIES
main.o		: main.cpp readlua.h

# HEADERS
HEADERS = readlua.h

readluafile: $(OBJS)
	${CC} -o $@ $^ $C $I $L

%.o: %.cpp $(HEADERS)
	${CC} -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o

