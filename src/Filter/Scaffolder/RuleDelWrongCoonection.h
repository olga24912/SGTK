#ifndef SCAFFOLDER_RULEDELWRONGCOONECTION_H
#define SCAFFOLDER_RULEDELWRONGCOONECTION_H

#include <set>
#include "Rule.h"

namespace filter {
    namespace scaffolder {
        class RuleDelWrongConnection : public Rule {
        private:
            std::string coordFile = "out.coord";
        public:
            void simplifyGraph(ContigGraph *filter) override;

        };
    }
}

#endif //SCAFFOLDER_RULEDELWRONGCOONECTION_H
