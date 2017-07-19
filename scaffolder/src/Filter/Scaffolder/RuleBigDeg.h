#ifndef SCAFFOLDER_RULEBIGDEG_H
#define SCAFFOLDER_RULEBIGDEG_H

#include <ContigGraph/ContigGraph.h>
#include "Rule.h"

namespace filter {
    namespace scaffolder {
        using namespace contig_graph;
        class RuleBigDeg : public Rule {
        private:
            const int BIG_DEG = 20;
            const int MAX_WEIGHT = 4;
        public:
            void simplifyGraph(ContigGraph *graph) override;

            void ignoreEdges(ContigGraph *graph, std::vector<int> vector);

            void ignoreEdgesR(ContigGraph *graph, std::vector<int> vector);
        };
    }
}


#endif //SCAFFOLDER_RULEBIGDEG_H
