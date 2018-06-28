#include "GFAGraphBuilder.h"

namespace builder {
    namespace graph_builder {
        void GFAGraphBuilder::setGFAFile(std::string file_name) {
            gfaFileName = file_name;

        }

        void GFAGraphBuilder::evaluate() {
            parseGFA();
        }

        void GFAGraphBuilder::handleS(std::string& name, std::string& record, bool firstLib, int i) {
            contigsId[name] = 2 * i;
            contigsName.push_back(name);

            if (firstLib) {
                graph->addVertex(2 * i, name, (int)record.size());
            }

            name += "-rev";

            contigsId[name] = 2 * i + 1;
            contigsName.push_back(name);

            if (firstLib) {
                graph->addVertex(2 * i + 1, name, (int)record.size());
            }
        }

        void GFAGraphBuilder::parseGFA() {
            bool firstLib = (graph->getLibNum() == 1);
            std::ifstream fin(gfaFileName);
            std::string s;
            int i = 0;
            while (getline(fin, s)) {
                std::stringstream ss(s);
                char RecordType;
                ss >> RecordType;
                if (RecordType == 'S') {
                    std::string name;
                    ss >> name;
                    std::string record;
                    ss >> record;
                    handleS(name, record, firstLib, i);
                    ++i;
                } else if (RecordType == 'L') {
                    std::string n1, n2, o1, o2;
                    ss >> n1 >> o1 >> n2 >> o2;
                    if (o1 == "-") {
                        n1 += "-rev";
                    }
                    if (o2 == "-") {
                        n2 += "-rev";
                    }

                    graph->addEdge(contigsId[n1], contigsId[n2], 1, 0, "");
                    graph->addEdge(contigsId[n2]^1, contigsId[n1]^1, 1, 0, "");
                }
            }

            fin.close();
        }
    }
}