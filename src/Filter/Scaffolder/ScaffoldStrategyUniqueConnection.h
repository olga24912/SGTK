#ifndef SCAFFOLDER_SCAFFOLDSTRATEGYUNIQUECONNECTION_H
#define SCAFFOLDER_SCAFFOLDSTRATEGYUNIQUECONNECTION_H

#include <Filter/ContigGraph/ContigGraph.h>
#include "ScaffoldStrategy.h"

namespace filter {
    namespace scaffolder {
        using namespace contig_graph;
        class ScaffoldStrategyUniqueConnection : public ScaffoldStrategy {
        private:
            bool isUniquePair(int v1, int v2, ContigGraph *graph);

        public:
            void addConnection(Scaffolds *scaffolds, ContigGraph *graph) override;

        private:
            DECL_LOGGER("ScaffoldStrategyUniqueConnection");
        };
    }
}

#endif //SCAFFOLDER_SCAFFOLDSTRATEGYUNIQUECONNECTION_H
