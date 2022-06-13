#ifndef MOLECULE_H
#define MOLECULE_H

#include "Atom.hpp"
#include <string>

class Molecule {
public:
    std::string dir;
    int natom;
    Atom *atoms;
    
    Molecule(std::string dir);
    ~Molecule();
    void read_coord(); 
};

#endif /* MOLECULE_H */