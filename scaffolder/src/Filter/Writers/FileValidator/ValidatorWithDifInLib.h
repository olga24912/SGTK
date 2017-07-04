#ifndef SCAFFOLDER_VALIDATORWITHDIFINLIB_H
#define SCAFFOLDER_VALIDATORWITHDIFINLIB_H

#include "FileValidator.h"

namespace filter {
    namespace writers {
        class ValidatorWithDifInLib : public FileValidator {
        private:
            std::vector<int> libs;
        public:
            ValidatorWithDifInLib(std::vector<int> libs) : libs(libs) {}

            bool isGoodVertexSet(std::vector<int> vert, Filter *filter) override;
        };
    }
}

#endif //SCAFFOLDER_VALIDATORWITHDIFINLIB_H
