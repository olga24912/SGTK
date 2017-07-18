#ifndef SCAFFOLDER_RULEDELSMALLCYCLE_H
#define SCAFFOLDER_RULEDELSMALLCYCLE_H

#include "Rule.h"

namespace filter {
    namespace scaffolder {
        class RuleDelSmallCycle : public Rule {
        private:
            const int MAX_DIF = 5;
        public:
            void simplifyGraph(ContigGraph *filter) override;

        };
    }
}

#endif //SCAFFOLDER_RULEDELSMALLCYCLE_H
