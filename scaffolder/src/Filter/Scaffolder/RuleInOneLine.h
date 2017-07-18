#ifndef SCAFFOLDER_RULEINONELINE_H
#define SCAFFOLDER_RULEINONELINE_H

#include "Rule.h"

namespace filter {
    namespace scaffolder {
        class RuleInOneLine : Rule {
        public:
            void simplifyGraph(ContigGraph *filter) override;

        private:
            bool havePath(ContigGraph *filter, int u, int w);
        };
    }
}


#endif //SCAFFOLDER_RULEINONELINE_H
