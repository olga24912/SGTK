#include "CommandAddInfoToGraph.h"

namespace filter {
    namespace commands {
        using namespace contig_graph;
        void CommandAddInfoToGraph::execute(std::string argv, State &state, ContigGraph &graph) {
            std::stringstream ss(argv);
            std::string infoFileName;
            std::string libName;
            std::string color;

            ss >> infoFileName >> libName >> color;

            INFO("start add info file \""<< infoFileName << "\" with new lib name \"" << libName << "\" and with new lib color " << color);


            std::ifstream infoin(infoFileName);
            int lib = graph.addLib(color, libName, ContigGraph::Lib::Type::SCAFF);

            std::string cur;
            while (getline(infoin, cur)) {
                std::stringstream ss(cur);
                std::string x;
                ss >> x;
                std::string contigName;
                int id;
                char dir;
                char tr;

                int v = -1;
                while (ss >> contigName >> id >> dir >> tr) {
                    contigName = contigName.substr(1);
                    int u = graph.getTargetId(contigName);
                    if (dir == '-') {
                        u ^= 1;
                    }

                    if (v != -1) {
                        graph.addEdge(v, u, lib, 1, 0, 0, 0, 0);
                        graph.addEdge(u^1, v^1, lib, 1, 0, 0, 0, 0);
                    }

                    v = u;
                }
            }

            infoin.close();
        }
    }
}