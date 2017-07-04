#ifndef SCAFFOLDER_SCAFFOLDSTRATEGYUNIQUECONNECTION_H
#define SCAFFOLDER_SCAFFOLDSTRATEGYUNIQUECONNECTION_H

#include "ScaffoldStrategy.h"

namespace filter {
    namespace scaffolder {
        class ScaffoldStrategyUniqueConnection : public ScaffoldStrategy {
        private:
            bool isUniquePair(int v1, int v2, Filter *graph);

        public:
            void addConnection(Scaffolds *scaffolds, Filter *graph, std::vector<int> minW) override;

        private:
            DECL_LOGGER("ScaffoldStrategyUniqueConnection");
        };
    }
}

#endif //SCAFFOLDER_SCAFFOLDSTRATEGYUNIQUECONNECTION_H
