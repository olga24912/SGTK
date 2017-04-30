#ifndef SCAFFOLDER_VALIDATORFORK_H
#define SCAFFOLDER_VALIDATORFORK_H


#include "FileValidator.h"

class ValidatorFork : public FileValidator {
private:
    int lib;
public:
    ValidatorFork(int lib) : lib(lib) {}

    bool isGoodVertexSet(std::vector<int> vert, Filter *filter);
};


#endif //SCAFFOLDER_VALIDATORFORK_H
