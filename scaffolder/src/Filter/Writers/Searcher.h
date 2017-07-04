#ifndef SCAFFOLDER_SEARCHER_H
#define SCAFFOLDER_SEARCHER_H

#include "Filter/Filters/Filter.h"
#include <algorithm>

namespace filter {
    namespace writers {
//class for search in graph.
        class Searcher {
        private:
            Filter *filter;

            void dfsFindComponent(int v, int currentCol, int *color);

        public:
            Searcher(Filter *filter1) : filter(filter1) {}

            //find all vertex on distance less then dist from v
            std::vector<int> findVertInLocalArea(int v, int dist);

            //return number of commponent, col[v] - #commponent for vertex v
            int findComponent(int *col);

        private:
            DECL_LOGGER("Searcher");
        };
    }
}

#endif //SCAFFOLDER_SEARCHER_H
