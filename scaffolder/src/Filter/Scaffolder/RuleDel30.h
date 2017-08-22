#ifndef SCAFFOLDER_RULEDELUNSUPPORTED30_H
#define SCAFFOLDER_RULEDELUNSUPPORTED30_H

#include "Rule.h"

namespace filter {
    namespace scaffolder {
        class RuleDel30 : Rule {
        public:
            void simplifyGraph(ContigGraph *filter) override;
        };
    }
}


#endif //SCAFFOLDER_RULEDELUNSUPPORTED30_H
