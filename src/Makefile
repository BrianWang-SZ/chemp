chemp: main.o Molecule.o GeomSolver.o VibSolver.o HfSolver.o CCSolver.o EnergySolver.o Helper.o
	g++ main.o Molecule.o GeomSolver.o VibSolver.o HfSolver.o CCSolver.o EnergySolver.o Helper.o -o chemp

Molecule.o: Molecule.cpp
	g++ -c Molecule.cpp

GeomSolver.o: GeomSolver.cpp
	g++ -c GeomSolver.cpp

VibSolver.o: VibSolver.cpp
	g++ -c VibSolver.cpp

HfSolver.o: HfSolver.cpp
	g++ -c HfSolver.cpp

CCSolver.o: CCSolver.cpp
	g++ -c CCSolver.cpp

EnergySolver.o: EnergySolver.cpp
	g++ -c EnergySolver.cpp

Helper.o: Helper.cpp
	g++ -c Helper.cpp

main.o:
	g++ -c main.cpp

clean:
	rm *.o chemp
