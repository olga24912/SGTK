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
            if (this->contigFileName != "*" &&
                    this->contigFileName != "") {
                initGraph();
            }

        }
    }
}
