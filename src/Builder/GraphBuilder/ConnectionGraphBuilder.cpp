#include "ConnectionGraphBuilder.h"

namespace builder {
    namespace graph_builder {
        void ConnectionGraphBuilder::setConnectionFile(std::string file_name) {
            connectionFileName = file_name;
        }

        void ConnectionGraphBuilder::evaluate() {
            initGraph();
            parseConnection();
        }

        void ConnectionGraphBuilder::parseConnection() {
            std::ifstream in(connectionFileName);

            std::string fr, frd, sc, scd;
            double w;
            int len;

            while (in >> fr >> frd >> sc >> scd >> w >> len) {
                std::string info;
                getline(in, info);
                if (info != "") {
                    int i = 0;
                    while (i < info.size() && info[i] == ' ') {
                        ++i;
                    }
                    info = info.substr(i, info.size() - i);
                    while (info.size() > 0 && info[info.size() - 1] == ' ') {
                        info.resize(info.size() - 1);
                    }
                    if (info[0] == '"' && info[info.size() - 1] == '"') {
                        info = info.substr(1, info.size() - 2);
                    }
                }

                fr = fr.substr(1, fr.size() - 1);
                sc = sc.substr(1, sc.size() - 1);

                if (frd[0] == '-') {
                    fr += "-rev";
                }
                if (scd[0] == '-') {
                    sc += "-rev";
                }

                graph->addEdge(graph->getTargetId(fr), graph->getTargetId(sc), w, len, info);
                graph->addEdge(graph->getTargetId(sc)^1, graph->getTargetId(fr)^1, w, len, info);
            }

            in.close();
        }
    }
}
