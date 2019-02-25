//
// Created by olga on 25.02.19.
//

#ifndef SGTK_GFA2GRAPHBUILDER_H
#define SGTK_GFA2GRAPHBUILDER_H

#include "GraphBuilder.h"

namespace builder {
    namespace graph_builder {
        class GFA2GraphBuilder : public GraphBuilder {
        private:
            std::string gfaFileName;
        public:
            void setGFAFile(std::string file_name);
            void evaluate() override;
        private:
            DECL_LOGGER("GFA2GraphBuilder");
        };
    }
}


#endif //SGTK_GFA2GRAPHBUILDER_H
