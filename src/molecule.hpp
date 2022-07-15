#ifndef MOLECULE_H
#define MOLECULE_H

#include "Atom.hpp"
#include <string>

class Molecule {
public:
    std::string dir;
    int natom, nomo;
    Atom *atoms;
    
    Molecule(std::string dir);
    ~Molecule();
    void read_coord();
    int calc_nomo();
};

#endif /* MOLECULE_H */