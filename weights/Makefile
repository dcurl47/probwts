SHELL      = /bin/sh

CPP=g++
CC=g++
LINK=g++

OPTS= -O3 -funroll-loops -g
#OPTS= -g

LIBRARIES=-lm

EXECNAME="weights.x"

# objects to compile
OBJ = calcWeights.o 

# General rule for compiling
all: usual

clean: 
	rm *.o $(EXECNAME)

usual: $(OBJ)
	$(LINK) $(OPTS) -o $(EXECNAME) $(OBJ) $(LIBRARIES)

#list of rules for each individual file
calcWeights.o: calcWeights.cpp
	$(CPP) $(OPTS) -c -o calcWeights.o calcWeights.cpp

