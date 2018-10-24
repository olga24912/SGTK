#ifndef SCAFFOLDER_COMMANDDELCOORD_H
#define SCAFFOLDER_COMMANDDELCOORD_H

#include "Command.h"

namespace filter {
    namespace commands {
        /*
         * delCoord <lib>
         *
         * Merge edges between same nodes, weight = sum of weights, delete coordinate.
         */
        class CommandDelCoord : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override {
                int lib;
                std::stringstream ss(argv);
                ss >> lib;
                std::vector<int> vers = graph.getVertexList();
                for (int i = 0; i < (int)vers.size(); ++i) {
                    int v = vers[i];
                    std::vector<int> edges = graph.getEdges(v);
                    std::vector<int> del(edges.size(), 0);
                    for (int j = 0; j < (int)edges.size(); ++j) {
                        if (del[j]) continue;
                        if (graph.getEdgeLib(edges[j]) != lib) continue;
                        graph.setCoord(edges[j], 0, 0, 0, 0);
                        for (int g = j + 1; g < (int)edges.size(); ++g) {
                            if (graph.getEdgeTo(edges[j]) == graph.getEdgeTo(edges[g])) {
                                graph.setWeight(edges[j], graph.getEdgeWeight(edges[j]) + graph.getEdgeWeight(edges[g]));
                                graph.delEdge(edges[g]);
                                del[g] = 1;
                            }
                        }
                    }
                }

            }
        };
    }
}

#endif //SCAFFOLDER_COMMANDDELCOORD_H
