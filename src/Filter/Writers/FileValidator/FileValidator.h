#ifndef SCAFFOLDER_FILEVALIDATOR_H
#define SCAFFOLDER_FILEVALIDATOR_H

#include <vector>
#include "Filter/ContigGraph/ContigGraph.h"

namespace filter {
    namespace writers {
        using namespace contig_graph;
        class FileValidator {
        public:
            virtual bool isGoodVertexSet(std::vector<int> vert, ContigGraph *graph) {
                return true;
            }

        protected:
            DECL_LOGGER("FileValidator");
        };
    }
}

#endif //SCAFFOLDER_FILEVALIDATOR_H
