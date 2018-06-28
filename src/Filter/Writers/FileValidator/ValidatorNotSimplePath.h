#ifndef SCAFFOLDER_VALIDATORNOTSIMPLEPATH_H
#define SCAFFOLDER_VALIDATORNOTSIMPLEPATH_H

#include <Filter/ContigGraph/ContigGraph.h>
#include "FileValidator.h"

namespace filter {
    namespace writers {
        using namespace contig_graph;

        class ValidatorNotSimplePath : public FileValidator {
        public:
            bool isGoodVertexSet(std::vector<int> vert, ContigGraph *graph);
        };
    }
}


#endif //SCAFFOLDER_VALIDATORNOTSIMPLEPATH_H
