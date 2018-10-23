#ifndef SCAFFOLDER_RNASPLITGRAPHBUILDER_H
#define SCAFFOLDER_RNASPLITGRAPHBUILDER_H

#include "PairReadGraphBuilder.h"

namespace builder {
    namespace graph_builder {
        /*
         * Generate connection between contigs by RNA split reads alignments
         */
        class RNASplitGraphBuilder : public PairReadGraphBuilder {
        protected:
            int get2Target(const seqan::BamAlignmentRecord &read) const override;

            void incEdgeWeight(seqan::BamAlignmentRecord& read1, seqan::BamAlignmentRecord& read2) override;

            int get1Target(const seqan::BamAlignmentRecord &read) const override;

        private:
            std::pair<int, int> getCoord(seqan::BamAlignmentRecord& read, int target);

            int changeEdges(int v1, std::pair<int, int> c1, int v2, std::pair<int, int> c2);

            bool isGoodEdgeFor1(ContigGraph::Edge& edge, std::pair<int, int> c) override;

            bool isGoodEdgeFor2(ContigGraph::Edge& edge, std::pair<int, int> c) override ;

            std::pair<int, int> relaxCoord(std::pair<int, int> c1, std::pair<int, int> c);

            DECL_LOGGER("RNASplitGraphBuilder");
        };
    }
}

#endif //SCAFFOLDER_RNASPLITGRAPHBUILDER_H
