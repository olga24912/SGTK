#ifndef SCAFFOLDER_GRAPHSPLITTER_H
#define SCAFFOLDER_GRAPHSPLITTER_H

#include <vector>
#include "Filter/ContigGraph/ContigGraph.h"

namespace filter {
    namespace writers {
        using namespace contig_graph;
//split graph on small parts
        class GraphSplitter {
        private:
            int maxEdge = 40; //max edges in one component
            int maxVert = 20; //max vertexs in one component

            //0 - vertex was not used, 1 - current component, 2 - vertex was used before, 3 - was not used,
            //but it is bed vertex for this component
            std::vector<int> used;
            std::vector<std::vector<int> > res; //res[i] - vertexs in i's component.
            std::vector<int> edgeCol; //component for i-s edge.

            ContigGraph *graph;
            std::vector<int> vert;

            void clear();

            void findNewComp(int v);

        public:
            GraphSplitter() {}

            GraphSplitter(int maxVert, int maxEdge) : maxVert(maxVert), maxEdge(maxEdge) {}

            //in graph graph split this set of vertex
            std::vector<std::vector<int> > split(ContigGraph *graph1, std::vector<int> vert);

        private:
            DECL_LOGGER("GraphSplitter");
        };
    }
}
#endif //SCAFFOLDER_GRAPHSPLITTER_H
