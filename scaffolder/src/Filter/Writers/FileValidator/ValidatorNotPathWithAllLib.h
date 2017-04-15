#ifndef SCAFFOLDER_VALIDATORNOTPATHWITHALLLIB_H
#define SCAFFOLDER_VALIDATORNOTPATHWITHALLLIB_H

#include "FileValidator.h"

class ValidatorNotPathWithAllLib : public FileValidator {
public:
    virtual bool isGoodVertexSet(std::vector<int> vert, Filter* filter);
};


#endif //SCAFFOLDER_VALIDATORNOTPATHWITHALLLIB_H
