#include <Filter/CommandParsers/State.h>
#include "CommandMinEdgeWeight.h"

namespace filter {
    namespace commands {
        void CommandMinEdgeWeight::execute(std::string argv, State &state, ContigGraph &graph) {
            INFO("set ContigGraph min Edge weight");
            int libNum, w;
            std::stringstream ss(argv);
            ss >> libNum >> w;

            std::vector<int> vert = graph.getVertexList();
            for (int v : vert) {
                std::vector<int> edges = graph.getEdges(v);
                for (int e : edges) {
                    if (graph.getEdgeLib(e) == libNum && graph.getEdgeWeight(e) < w) {
                        graph.delEdge(e);
                    }
                }
            }
        }
    }
}
