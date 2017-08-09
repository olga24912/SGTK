#ifndef SCAFFOLDER_RNASPLITGRAPHBUILDER_H
#define SCAFFOLDER_RNASPLITGRAPHBUILDER_H

#include "RNAPairReadGraphBuilder.h"

namespace builder {
    namespace graph_builder {
        class RNASplitGraphBuilder : public RNAPairReadGraphBuilder {
        protected:
            int DIRECT = 0;

            int get2Target(const seqan::BamAlignmentRecord &read) const override;

            void incEdgeWeight(seqan::BamAlignmentRecord read1, seqan::BamAlignmentRecord read2) override;

            int get1Target(const seqan::BamAlignmentRecord &read) const override;

        private:
            DECL_LOGGER("RNASplitGraphBuilder");

            std::pair<int, int> getCoord(seqan::BamAlignmentRecord read);

            int changeEdges(int v1, std::pair<int, int> c1, int v2, std::pair<int, int> c2);

            bool isGoodEdgeFor1(ContigGraph::Edge edge, std::pair<int, int> c);

            bool isGoodEdgeFor2(ContigGraph::Edge edge, std::pair<int, int> c);

            std::pair<int, int> relaxCoord(std::pair<int, int> c1, std::pair<int, int> c);

        public:
            const int BACK = 1;
            const int FORWARD = 0;

            void setDirect(int dir) {
                DIRECT = dir;
            };
        };
    }
}

#endif //SCAFFOLDER_RNASPLITGRAPHBUILDER_H
