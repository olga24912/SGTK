#ifndef SCAFFOLDER_RNASPLIT30GRAPHBUILDER_H
#define SCAFFOLDER_RNASPLIT30GRAPHBUILDER_H

#include "RNASplitGraphBuilder.h"

namespace builder {
    namespace graph_builder {
        class RNASplit30GraphBuilder : public RNASplitGraphBuilder {
        protected:
            ContigGraph::Lib::Type getLibType() override;

        private:
            std::string getLibColor() override;
        private:
            DECL_LOGGER("RNASplit30GraphBuilder");
        };
    }
}


#endif //SCAFFOLDER_RNASPLIT30GRAPHBUILDER_H
