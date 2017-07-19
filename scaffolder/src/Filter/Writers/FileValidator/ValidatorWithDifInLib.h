#ifndef SCAFFOLDER_VALIDATORWITHDIFINLIB_H
#define SCAFFOLDER_VALIDATORWITHDIFINLIB_H

#include <ContigGraph/ContigGraph.h>
#include "FileValidator.h"

namespace filter {
    namespace writers {
        using namespace contig_graph;
        class ValidatorWithDifInLib : public FileValidator {
        private:
            std::vector<int> libs;
        public:
            ValidatorWithDifInLib(std::vector<int> libs) : libs(libs) {}

            bool isGoodVertexSet(std::vector<int> vert, ContigGraph *graph) override;
        };
    }
}

#endif //SCAFFOLDER_VALIDATORWITHDIFINLIB_H
