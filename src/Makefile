chemp: main.o Molecule.o GeomSolver.o VibSolver.o
	g++ main.o Molecule.o GeomSolver.o VibSolver.o -o chemp

Molecule.o: Molecule.cpp
	g++ -c Molecule.cpp

GeomSolver.o: GeomSolver.cpp
	g++ -c GeomSolver.cpp

VibSolver.o: VibSolver.cpp
	g++ -c VibSolver.cpp

main.o:
	g++ -c main.cpp

clean:
	rm *.o chemp
