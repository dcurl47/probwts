SHELL      = /bin/sh

CPP=g++
CC=g++
LINK=g++

OPTS= -O3 -funroll-loops -g
#OPTS= -g

LIBRARIES=-lm

EXECNAME="nne.x"

# objects to compile
OBJ = calcnne.o 

# General rule for compiling
all: usual

clean: 
	rm *.o $(EXECNAME)

usual: $(OBJ)
	$(LINK) $(OPTS) -o $(EXECNAME) $(OBJ) $(LIBRARIES)

#list of rules for each individual file
calcnne.o: calcnne.cpp
	$(CPP) $(OPTS) -c -o calcnne.o calcnne.cpp

