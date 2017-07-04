#ifndef SCAFFOLDER_DNAPAIRREADGRAPHBUILDER_H
#define SCAFFOLDER_DNAPAIRREADGRAPHBUILDER_H

#include "PairReadGraphBuilder.h"
namespace builder {
    namespace graph_builder {
//find conection between contigs by DNA pair reads
        class DNAPairReadGraphBuilder : public PairReadGraphBuilder {
        private:
            int distBetweenPairReads = (int) 1e9;
            std::map<std::string, int> read1DistToEnd;
            std::map<std::string, int> read2DistToEnd;

            void addInfoAboutRead(std::string readName, int target, seqan::BamAlignmentRecord read);

            void addInfoAbout2Read(std::string readName, int target, seqan::BamAlignmentRecord read);

            void incEdgeWeight(seqan::BamAlignmentRecord read1, seqan::BamAlignmentRecord read2);

            int readDist(seqan::BamAlignmentRecord read);

            virtual std::string getLibColor();

        protected:
            ContigGraph::Lib::Type getLibType() override;

        public:
            // set max distance between DNA pair read,
            // if  dist will be more this conatcion will be ignore
            void setDistBetweenPairReads(int distBetweenPairReads);

        private:
            DECL_LOGGER("DNAPairReadGraphBuilder");
        };
    }
}
#endif //SCAFFOLDER_DNAPAIRREADGRAPHBUILDER_H
