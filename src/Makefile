CC=g++
CFLAGS=-g -c -Wall
LDFLAGS=-g
SOURCES=main.cpp HFSolver.cpp GeomSolver.cpp VibSolver.cpp Mp2Solver.cpp Molecule.cpp Helper.cpp DIIS.cpp EnergySolver.cpp CCSolver.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=chemp

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm $(OBJECTS) chemp
