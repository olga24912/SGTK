//
// Created by olga on 25.02.19.
//

#include "GFA2GraphBuilder.h"

namespace builder {
    namespace graph_builder {
        void GFA2GraphBuilder::setGFAFile(std::string file_name) {
            this->gfaFileName = file_name;
        }

        void GFA2GraphBuilder::evaluate() {
            bool firstLib = (graph->getLibNum() == 1);
            if (this->contigFileName != "*" &&
                this->contigFileName != "" &&
                    !firstLib) {
                initGraph();
                firstLib = false;
            }

            std::ifstream fin(gfaFileName);
            std::string s;

            int i = 0;
            while (getline(fin, s)) {
                std::stringstream ss(s);
                char RecordType;
                ss >> RecordType;
                if (RecordType == 'S') {
                    if (!firstLib) {
                        std::string name;
                        ss >> name;
                        int slen;
                        ss >> slen;
                        std::string record;
                        ss >> record;

                        graph->addVertex(2 * i, name, slen);
                        graph->addVertex(2 * i + 1, name + "-rev", slen);
                    }
                    ++i;
                } else if (RecordType == 'E') {
                    std::string eid, sid1, sid2, beg1, end1, beg2, end2;
                    ss >> eid >> sid1 >> sid2 >> beg1 >> end1 >> beg2 >> end2;
                    if (sid1.back() == '-') {
                        sid1 += "rev";
                    } else {
                        sid1.pop_back();
                    }

                    if (sid2.back() == '-') {
                        sid2 += "rev";
                    } else {
                        sid2.pop_back();
                    }

                    if (end1.back() == '$' && beg2 == "0") {
                        graph->addEdge(graph->getTargetId(sid1), graph->getTargetId(sid2), 1, 0, "");
                        graph->addEdge(graph->getTargetId(sid2)^1, graph->getTargetId(sid1)^1, 1, 0, "");
                    } else if (end2.back() == '$' && beg1 == "0") {
                        graph->addEdge(graph->getTargetId(sid2), graph->getTargetId(sid1), 1, 0, "");
                        graph->addEdge(graph->getTargetId(sid1)^1, graph->getTargetId(sid2)^1, 1, 0, "");
                    }
                } else if (RecordType == 'G') {
                    std::string gid, sid1, sid2;
                    int dist;
                    ss >> gid >> sid1 >> sid2 >> dist;
                    if (sid1.back() == '-') {
                        sid1 += "rev";
                    } else {
                        sid1.pop_back();
                    }

                    if (sid2.back() == '-') {
                        sid2 += "rev";
                    } else {
                        sid2.pop_back();
                    }

                    graph->addEdge(graph->getTargetId(sid1), graph->getTargetId(sid2), 1, dist, "");
                    graph->addEdge(graph->getTargetId(sid2)^1, graph->getTargetId(sid1)^1, 1, dist, "");
                }
            }
        }
    }
}
