SHELL      = /bin/sh

CPP=g++
CC=g++
LINK=g++

OPTS= -O3 -funroll-loops -g
#OPTS= -g

LIBRARIES=-lm

EXECNAME="poly.x"

# objects to compile
OBJ = calcpoly.o 

# General rule for compiling
all: usual

clean: 
	rm *.o $(EXECNAME)

usual: $(OBJ)
	$(LINK) $(OPTS) -o $(EXECNAME) $(OBJ) $(LIBRARIES)

#list of rules for each individual file
calcpoly.o: calcpoly.cpp
	$(CPP) $(OPTS) -c -o calcpoly.o calcpoly.cpp

