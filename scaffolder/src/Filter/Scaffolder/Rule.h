#ifndef SCAFFOLDER_RULE_H
#define SCAFFOLDER_RULE_H

#include <ContigGraph/ContigGraph.h>

namespace filter {
    namespace scaffolder {
        using namespace contig_graph;
        class Rule {
        public:
            virtual void simplifyGraph(ContigGraph* filter) = 0;
        };
    }
}


#endif //SCAFFOLDER_RULE_H
