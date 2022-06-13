#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "Molecule.hpp"

Molecule::Molecule(std::string dir){
    this -> dir = dir;
    read_coord();
}

Molecule::~Molecule(){
    delete[] atoms;
}

void Molecule::read_coord(){
    std::string path = dir + "/geom.dat";
    // open input file
    FILE *in;
    
    if ((in = fopen(path.c_str(), "r")) == NULL) {
        perror("Error opening geom file");
        exit(-1);
    }
    
    // read in the number of atoms
    int num;
    if (fscanf(in, "%d", &num) != 1) {       
        perror("Error reading number of atoms");
        exit(-1);
    }
    natom = num;

    // read each atom
    int count = 0;
    
    Atom *atoms = new Atom[natom];

    double zval = 0.0, x = 0.0,
              y = 0.0, z = 0.0;
    
    while (fscanf(in, "%lf %lf %lf %lf", &zval, &x, &y, &z) == 4) {
        
        atoms[count].zval = (int) zval;
        atoms[count].x = x;
        atoms[count].y = y;
        atoms[count].z = z;

        count++;
    }

    if (count < num) {
        perror("Atom number does not match data");
        exit(-1);
    }

    if (fclose(in)){
        perror("Error closing file stream!");
        exit(-1);
    }

    this -> atoms = atoms;
}
