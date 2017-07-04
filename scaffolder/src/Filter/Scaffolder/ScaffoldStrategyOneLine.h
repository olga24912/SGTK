#ifndef SCAFFOLDER_SCAFFOLDSTRATEGYONELINE_H
#define SCAFFOLDER_SCAFFOLDSTRATEGYONELINE_H

#include "ScaffoldStrategy.h"

namespace filter {
    namespace scaffolder {
        class ScaffoldStrategyOneLine : public ScaffoldStrategy {
        protected:
            std::vector<int> newCon;
            std::vector<int> newConR;

            void addFirstConnection(Scaffolds *scaffolds, Filter *graph, std::vector<int> minW);

            void delEdgeFromDifPath(Scaffolds *scaffolds, Filter *graph);

        public:
            void addConnection(Scaffolds *scaffolds, Filter *graph, std::vector<int> minW) override;

        private:
            DECL_LOGGER("ScaffoldStrategyOneLine");
        };
    }
}

#endif //SCAFFOLDER_SCAFFOLDSTRATEGYONELINE_H
