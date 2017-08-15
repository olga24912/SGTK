#ifndef SCAFFOLDER_RULE_H
#define SCAFFOLDER_RULE_H

#include <ContigGraph/ContigGraph.h>

namespace filter {
    namespace scaffolder {
        using namespace contig_graph;
        class Rule {
        public:
            virtual void simplifyGraph(ContigGraph* filter) = 0;

            bool sameCoord1(ContigGraph *graph, int e1, int e2);

            bool sameCoord2(ContigGraph *graph, int e1, int e2);
        };
    }
}


#endif //SCAFFOLDER_RULE_H
