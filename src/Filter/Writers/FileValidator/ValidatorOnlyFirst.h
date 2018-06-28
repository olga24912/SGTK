#ifndef SCAFFOLDER_VALIDATORONLYFIRST_H
#define SCAFFOLDER_VALIDATORONLYFIRST_H

#include <vector>
#include <ContigGraph/ContigGraph.h>
#include "FileValidator.h"

namespace filter {
    namespace writers {
        using namespace contig_graph;
        class ValidatorOnlyFirst : public FileValidator {
        private:
            int libPr, libNotPr;
        public:
            ValidatorOnlyFirst(int libPr, int libNotPr) : libPr(libPr), libNotPr(libNotPr) {}

            bool isGoodVertexSet(std::vector<int> vert, ContigGraph *graph) override;
        };
    }
}

#endif //SCAFFOLDER_VALIDATORONLYFIRST_H
