#ifndef SCAFFOLDER_DNAMATEPAIRREADSGRAPHBUILDER_H
#define SCAFFOLDER_DNAMATEPAIRREADSGRAPHBUILDER_H


#include "PairReadGraphBuilder.h"

namespace builder {
    namespace graph_builder {
        class DNAMatePairReadsGraphBuilder : public PairReadGraphBuilder {
        protected:
            int get2Target(const seqan::BamAlignmentRecord &read) const override;

            void incEdgeWeight(seqan::BamAlignmentRecord& read1, seqan::BamAlignmentRecord& read2) override;

            int get1Target(const seqan::BamAlignmentRecord &read) const override;

        private:
            DECL_LOGGER("DNAMatePairReadsGraphBuilder");
            int changeEdges(int v1, std::pair<int, int> c1, int v2, std::pair<int, int> c2);
        };
    }
}


#endif //SCAFFOLDER_DNAMATEPAIRREADSGRAPHBUILDER_H
