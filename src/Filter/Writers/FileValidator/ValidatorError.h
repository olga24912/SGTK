#ifndef SCAFFOLDER_VALIDATORERROR_H
#define SCAFFOLDER_VALIDATORERROR_H

#include <Filter/AligInfo/InfoAboutContigsAlig.h>
#include <ContigGraph/ContigGraph.h>
#include "FileValidator.h"

namespace filter {
    namespace writers {
        using namespace alig_info;
        using namespace contig_graph;
        class ValidatorError : public FileValidator {
        private:
            InfoAboutContigsAlig aligInfo;
            std::string coordFileName;
            int libError;
        public:
            ValidatorError(std::string coordFile, int libError) : coordFileName(coordFile), libError(libError) {};
            bool isGoodVertexSet(std::vector<int> vert, ContigGraph *graph) override;
        };
    }
}

#endif //SCAFFOLDER_VALIDATORERROR_H
