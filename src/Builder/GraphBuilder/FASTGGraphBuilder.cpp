#include "FASTGGraphBuilder.h"

namespace builder {
    namespace graph_builder {
        void FASTGGraphBuilder::setFASTGFile(std::string file_name) {
            fastgFileName = file_name;
        }

        void FASTGGraphBuilder::evaluate() {
            initGraph();
            parseFASTG();
        }

        std::string FASTGGraphBuilder::parseFirstEdge(std::string s) {
            std::string res="";
            int i =0;
            while (i < s.size() && s[i] != ':' && s[i] != ';') {
                res += s[i];
                ++i;
            }
            if (res[res.size() - 1] == '\'') {
                res.resize(res.size() - 1);
                res += "-rev";
            }
            return res;
        }

        std::vector<std::string> FASTGGraphBuilder::parseFollowingEdges(std::string s) {
            int i = 0;
            while (i < s.size() && s[i] != ':') {
                ++i;
            }
            ++i;
            std::vector<std::string> res;
            while (i < s.size()) {
                std::string cur = "";
                while (i < s.size() && s[i] != ',' && s[i] != ';') {
                    cur += s[i];
                    ++i;
                }
                if (cur[cur.size() - 1] == '\'') {
                    cur.resize(cur.size() - 1);
                    cur += "-rev";
                }
                res.push_back(cur);
                ++i;
            }
            return res;
        }

        void FASTGGraphBuilder::parseFASTG() {
            std::ifstream ifs(fastgFileName);
            std::string name;
            while (std::getline(ifs, name)) {
                if (name.size() > 0 && name[0] == '>') {
                    name = name.substr(1);
                    std::vector<std::string> followContigs = parseFollowingEdges(name);
                    std::string sc = parseFirstEdge(name);

                    for (int i = 0; i < followContigs.size(); ++i) {
                        graph->addEdge(graph->getTargetId(followContigs[i]), graph->getTargetId(sc), 1, 0, "");
                    }
                }
            }
        }
    }
}
