#ifndef SCAFFOLDER_SCAFFOLDSTRATEGYONELINE_H
#define SCAFFOLDER_SCAFFOLDSTRATEGYONELINE_H


#include "ScaffoldStrategy.h"

namespace filter {
    namespace scaffolder {
        using namespace contig_graph;

        class ScaffoldStrategyOneLine : public ScaffoldStrategy {
        protected:
            void addFirstConnection(Scaffolds *scaffolds, ContigGraph *graph);

            void delEdgeFromDifPath(Scaffolds *scaffolds, ContigGraph *graph);

        public:
            void addConnection(Scaffolds *scaffolds, ContigGraph *graph) override;

        };
    }
}


#endif //SCAFFOLDER_SCAFFOLDSTRATEGYONELINE_H
