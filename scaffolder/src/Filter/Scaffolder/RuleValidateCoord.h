#ifndef SCAFFOLDER_RULEVALIDATECOORD_H
#define SCAFFOLDER_RULEVALIDATECOORD_H

#include <Filter/AligInfo/InfoAboutContigsAlig.h>
#include "Rule.h"

namespace filter {
    namespace scaffolder {
        class RuleValidateCoord : public Rule {
        public:
            void simplifyGraph(ContigGraph *filter) override;

            void validateEdge(ContigGraph *graph, int e);
        };
    }
}


#endif //SCAFFOLDER_RULEVALIDATECOORD_H
