#ifndef ENERGYSOLVER_H
#define ENERGYSOLVER_H

#include <string>
#include "type.h"

class EnergySolver(){
public:
    virtual double compute() = 0;
}

#endif /* ENERGYSOLVER_H */