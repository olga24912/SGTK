#ifndef SCAFFOLDER_RNAPAIRREADGRAPHBUILDER_H
#define SCAFFOLDER_RNAPAIRREADGRAPHBUILDER_H

#include "PairReadGraphBuilder.h"

namespace builder {
    namespace graph_builder {
//generate connection between contigs by RNA pair reads
        class RNAPairReadGraphBuilder : public PairReadGraphBuilder {
        private:
            virtual std::string getLibColor();

        protected:
            ContigGraph::Lib::Type getLibType() override;

        public:
            RNAPairReadGraphBuilder() = default;

        private:
            DECL_LOGGER("RNAPairReadGraphBuilder");
        };
    }
}


#endif //SCAFFOLDER_RNAPAIRREADGRAPHBUILDER_H
