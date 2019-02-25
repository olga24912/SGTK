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
                        handleSline(name, slen, i);
                    }
                    ++i;
                } else if (RecordType == 'E') {
                } else if (RecordType == 'G') {
                }
            }



        }

        void GFA2GraphBuilder::handleEline() {

        }

        void GFA2GraphBuilder::handleGline() {

        }

        void GFA2GraphBuilder::handleSline(std::string &name, int slen, int i) {
            graph->addVertex(2 * i, name, slen);
            graph->addVertex(2 * i + 1, name + "-rev", slen);
        }
    }
}
