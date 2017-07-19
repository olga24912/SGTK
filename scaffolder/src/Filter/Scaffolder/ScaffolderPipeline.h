#ifndef SCAFFOLDER_SCAFFOLDERPIPELINE_H
#define SCAFFOLDER_SCAFFOLDERPIPELINE_H

#include <string>
#include <ContigGraph/ContigGraph.h>
#include "Scaffolds.h"

namespace filter {
    namespace scaffolder {
        using namespace contig_graph;
        class ScaffolderPipeline {
        private:

        public:
            void evaluate(ContigGraph *graph, std::string contigFile, std::string out);

        private:
            DECL_LOGGER("ScaffolderPipeline");
        };
    }
}
#endif //SCAFFOLDER_SCAFFOLDERPIPELINE_H
