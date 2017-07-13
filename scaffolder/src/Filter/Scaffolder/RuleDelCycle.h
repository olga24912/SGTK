#ifndef SCAFFOLDER_RULEDELCYCLE_H
#define SCAFFOLDER_RULEDELCYCLE_H

#include "Rule.h"

namespace filter {
    namespace scaffolder {
        class RuleDelCycle : public Rule {
        private:
            const int DEL_DIF = 2;

            std::vector<int> topsort;
            std::vector<int> color;
            std::vector<int> topSortPos;

            void topSortDfs(int v, Filter *graph, std::vector<int> *used);

            void colorDfs(int v, int col, Filter *graph);

            void topSort(Filter *graph);

            int findCycle(Filter *graph);
        public:
            void simplifyGraph(Filter *filter) override;
        };
    }
}
#endif //SCAFFOLDER_RULEDELCYCLE_H
