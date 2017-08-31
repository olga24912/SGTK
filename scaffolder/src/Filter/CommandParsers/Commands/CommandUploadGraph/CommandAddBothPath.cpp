#include "CommandAddBothPath.h"

namespace filter {
    namespace commands {
        using namespace contig_graph;
        void CommandAddBothPath::execute(std::string argv, State &state, contig_graph::ContigGraph &graph) {
            std::stringstream ss(argv);
            std::string fileName;
            std::string libName;
            std::string color;

            ss >> fileName >> libName >> color;

            std::ifstream infoin(fileName);
            int lib = graph.addLib(color, libName, ContigGraph::Lib::SCAFF);

            std::string cur;
            while (getline(infoin, cur)) {
                std::vector<int> vid;
                std::string curName = "";
                int i = 0;
                while (i <= (int)cur.size()) {
                    if (i < cur.size() && cur[i] != '-' && cur[i] != '>') {
                        curName += cur[i];
                    } else if (cur[i] == '-' || i == (int)cur.size()){
                        if (curName[curName.size() - 1] != ')') {
                            if (curName[curName.size() - 2] == '/') {
                                curName.resize(curName.size() - 2);
                                int v = graph.getTargetId(curName);
                                vid.push_back(v^1);
                            } else {
                                vid.push_back(graph.getTargetId(curName));
                            }
                        }
                        curName = "";
                    }
                    ++i;
                }
                for (int j = 1; j < (int)vid.size(); ++j) {
                    graph.addEdge(vid[j - 1], vid[j], lib, 1, 0, 0, 0, 0);
                    graph.addEdge(vid[j] ^ 1, vid[j - 1] ^ 1, lib, 1, 0, 0, 0, 0);
                }
            }

            infoin.close();

        }
    }
}
