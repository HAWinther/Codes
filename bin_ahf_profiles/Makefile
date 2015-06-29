#####################################################
# Hans A. Winther (2015) (hans.a.winther@gmail.com)
#####################################################

SHELL := /bin/bash

# Set compiler
# CC = g++-mp-4.8

# Includes for LUA 5.x or above
L = -L/Users/winther/local/lua-5.2.3/install/lib -L/Users/winther/local/lib -llua -ldl
I = -I/Users/winther/local/lua-5.2.3/install/include -I/Users/winther/local/include -lm

OPTIONS = -O3

#####################################################
# Number of dimensions and working precision (FLOAT,DOUBLE,LONGDOUBLE)
#####################################################
NDIM=3
TYPE=DOUBLE

#####################################################
# Profile limits rvir * RVIRMINFAC < r < rvir * RVIRMAXFAC
#####################################################
RVIRMINFAC=0.01
OPTIONS += -DRVIRMINFAC=$(RVIRMINFAC)

#####################################################
# Bin AHF profiles
#####################################################
NCOLPROFILEFILE=24
OPTIONS += -DNCOLPROFILEFILE=$(NCOLPROFILEFILE)

#####################################################
# Set all the options
#####################################################
OPTIONS += -Wall -DNDIM=$(NDIM) -D$(TYPE)
#####################################################

C = $(OPTIONS) 
TARGETS := calc_halo_profiles
all: $(TARGETS)

# OBJECT FILES
OBJS = main.o

# DEPENDENCIES
main.o		: global.h io.h analysis.h Makefile

# HEADERS
HEADERS = global.h analysis.h io.h

calc_halo_profiles: $(OBJS)
	$(CC) -o $@ $^ $C $I $L $(SFML)

%.o: %.cpp $(HEADERS)
	$(CC) -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o

