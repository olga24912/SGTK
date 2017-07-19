#ifndef SCAFFOLDER_VALIDATORNOTPATHWITHALLLIB_H
#define SCAFFOLDER_VALIDATORNOTPATHWITHALLLIB_H

#include "Filter/ContigGraph/ContigGraph.h"
#include "FileValidator.h"

namespace filter {
    namespace writers {
        using namespace contig_graph;
        class ValidatorNotPathWithAllLib : public FileValidator {
        public:
            virtual bool isGoodVertexSet(std::vector<int> vert, ContigGraph *graph);
        };
    }
}

#endif //SCAFFOLDER_VALIDATORNOTPATHWITHALLLIB_H
