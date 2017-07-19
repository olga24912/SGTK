#ifndef SCAFFOLDER_SEARCHER_H
#define SCAFFOLDER_SEARCHER_H

#include <algorithm>
#include "Filter/ContigGraph/ContigGraph.h"

using namespace filter::contig_graph;

namespace filter {
    namespace writers {
//class for search in graph.
        class Searcher {
        private:
            ContigGraph *graph;

            void dfsFindComponent(int v, int currentCol, int *color);

        public:
            Searcher(ContigGraph *graph1) : graph(graph1) {}

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
