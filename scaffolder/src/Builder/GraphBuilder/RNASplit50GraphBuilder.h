#ifndef SCAFFOLDER_RNASPLIT50GRAPHBUILDER_H
#define SCAFFOLDER_RNASPLIT50GRAPHBUILDER_H

#include "RNASplitGraphBuilder.h"

namespace builder {
    namespace graph_builder {
        class RNASplit50GraphBuilder : public RNASplitGraphBuilder {
        protected:
            ContigGraph::Lib::Type getLibType() override;

        private:
            std::string getLibColor() override;
        private:
            DECL_LOGGER("RNASplit50GraphBuilder");

        };
    }
}


#endif //SCAFFOLDER_RNASPLIT50GRAPHBUILDER_H
