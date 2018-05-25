#ifndef SCAFFOLDER_DNAMATEPAIRREADSGRAPHBUILDER_H
#define SCAFFOLDER_DNAMATEPAIRREADSGRAPHBUILDER_H


#include "PairReadGraphBuilder.h"

namespace builder {
    namespace graph_builder {
        class DNAMatePairReadsGraphBuilder : public PairReadGraphBuilder {
        protected:
            int firstRev = 0;
            int secondRev = 1;
            int get2Target(const seqan::BamAlignmentRecord &read) const override;
            void incEdgeWeight(seqan::BamAlignmentRecord& read1, seqan::BamAlignmentRecord& read2) override;
            int get1Target(const seqan::BamAlignmentRecord &read) const override;

        private:
            DECL_LOGGER("DNAMatePairReadsGraphBuilder");
            int changeEdges(int v1, std::pair<int, int> c1, int v2, std::pair<int, int> c2);

        public:
            void setRevFirstFlag(int firstRev) {
                this->firstRev = firstRev;
            }

            void setRevSecondFlag(int secondRev) {
                this->secondRev = secondRev;
            }
        };
    }
}


#endif //SCAFFOLDER_DNAMATEPAIRREADSGRAPHBUILDER_H
