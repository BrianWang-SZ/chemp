#ifndef MOLECULE_H
#define MOLECULE_H

#include "Atom.hpp"

class Molecule {
public:
    int natom;
    Atom *atoms;
    
    Molecule(const char*);
    ~Molecule();
    void read_coord(const char*);   
};

#endif /* MOLECULE_H */