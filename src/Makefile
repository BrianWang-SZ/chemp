chemp: main.o Molecule.o GeomSolver.o
	g++ main.o Molecule.o GeomSolver.o -o chemp

Molecule.o: Molecule.cpp
	g++ -c Molecule.cpp

GeomSolver.o: GeomSolver.cpp
	g++ -c GeomSolver.cpp

main.o:
	g++ -c main.cpp

clean:
	rm *.o chemp
