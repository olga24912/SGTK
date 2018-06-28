#ifndef SCAFFOLDER_VALIDATORFEWPARTS_H
#define SCAFFOLDER_VALIDATORFEWPARTS_H

#include <ContigGraph/ContigGraph.h>
#include "FileValidator.h"

namespace filter {
    namespace writers {
        using namespace contig_graph;
        class ValidatorFewParts : public FileValidator {
        public:
            bool isGoodVertexSet(std::vector<int> vert, ContigGraph *graph) override;
        };
    }
}

#endif //SCAFFOLDER_VALIDATORFEWPARTS_H
