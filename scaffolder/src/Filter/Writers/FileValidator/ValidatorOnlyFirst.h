#ifndef SCAFFOLDER_VALIDATORONLYFIRST_H
#define SCAFFOLDER_VALIDATORONLYFIRST_H

#include <vector>
#include "FileValidator.h"

namespace filter {
    namespace writers {
        class ValidatorOnlyFirst : public FileValidator {
        private:
            int libPr, libNotPr;
        public:
            ValidatorOnlyFirst(int libPr, int libNotPr) : libPr(libPr), libNotPr(libNotPr) {}

            bool isGoodVertexSet(std::vector<int> vert, Filter *filter) override;
        };
    }
}

#endif //SCAFFOLDER_VALIDATORONLYFIRST_H
