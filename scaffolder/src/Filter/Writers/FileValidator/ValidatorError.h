#ifndef SCAFFOLDER_VALIDATORERROR_H
#define SCAFFOLDER_VALIDATORERROR_H

#include <Filter/Statistics/InfoAboutContigsAlig.h>
#include "FileValidator.h"

namespace filter {
    namespace writers {
        using namespace statistics;
        class ValidatorError : public FileValidator {
        private:
            InfoAboutContigsAlig aligInfo;
            std::string coordFileName;
            int libError;
        public:
            ValidatorError(std::string coordFile, int libError) : coordFileName(coordFile), libError(libError) {};
            bool isGoodVertexSet(std::vector<int> vert, Filter *filter) override;
        };
    }
}

#endif //SCAFFOLDER_VALIDATORERROR_H
