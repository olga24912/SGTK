#ifndef SCAFFOLDER_RULE_H
#define SCAFFOLDER_RULE_H

#include <Filter/Filters/ContigGraph.h>

namespace filter {
    namespace scaffolder {
        class Rule {
        public:
            virtual void simplifyGraph(ContigGraph* filter) = 0;
        };
    }
}


#endif //SCAFFOLDER_RULE_H
