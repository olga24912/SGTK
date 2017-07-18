#ifndef SCAFFOLDER_VALIDATORFEWPARTS_H
#define SCAFFOLDER_VALIDATORFEWPARTS_H

#include "FileValidator.h"

namespace filter {
    namespace writers {
        class ValidatorFewParts : public FileValidator {
        public:
            bool isGoodVertexSet(std::vector<int> vert, ContigGraph *graph) override;
        };
    }
}

#endif //SCAFFOLDER_VALIDATORFEWPARTS_H
