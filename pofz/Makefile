SHELL      = /bin/sh

CPP=g++
CC=g++
LINK=g++

OPTS= -O3 -funroll-loops -g
#OPTS= -g

LIBRARIES=-lm

EXECNAME="pofz.x"

# objects to compile
OBJ = calcpofz.o 

# General rule for compiling
all: usual

clean: 
	rm *.o $(EXECNAME)

usual: $(OBJ)
	$(LINK) $(OPTS) -o $(EXECNAME) $(OBJ) $(LIBRARIES)

#list of rules for each individual file
calcpofz.o: calcpofz.cpp
	$(CPP) $(OPTS) -c -o calcpofz.o calcpofz.cpp

