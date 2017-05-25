#ifndef SCAFFOLDER_VALIDATORFEWPARTS_H
#define SCAFFOLDER_VALIDATORFEWPARTS_H


#include "FileValidator.h"

class ValidatorFewParts : public FileValidator {
public:
    bool isGoodVertexSet(std::vector<int> vert, Filter *filter) override;
};


#endif //SCAFFOLDER_VALIDATORFEWPARTS_H
