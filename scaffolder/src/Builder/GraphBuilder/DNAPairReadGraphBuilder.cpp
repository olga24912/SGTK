#include "DNAPairReadGraphBuilder.h"
#include "ReadsSplitter/Utils/SeqanUtils.h"

namespace builder {
    namespace graph_builder {
        void DNAPairReadGraphBuilder::setDistBetweenPairReads(int distBetweenPairReads) {
            TRACE("setDistBetweenPairReads distBetweenPairReads=" << distBetweenPairReads);
            DNAPairReadGraphBuilder::distBetweenPairReads = distBetweenPairReads;
        }

        /*void
        DNAPairReadGraphBuilder::addInfoAboutRead(std::string readName, int target, seqan::BamAlignmentRecord read) {
            PairReadGraphBuilder::addInfoAboutRead(readName, target, read);
            TRACE("addInfoAboutRead readName=" << readName << " v=" << target);
            read1DistToEnd[readName] = readDist(read);
        }*/

        int DNAPairReadGraphBuilder::readDist(seqan::BamAlignmentRecord read) {
            TRACE("get readDist to the end");
            if (!seqan::hasFlagRC(read)) {
                return (graph->getTargetLen(2 * read.rID) - read.beginPos -
                        (int) (read.seq.data_end - read.seq.data_begin));
            } else {
                return (read.beginPos);
            }
        }

        /*void
        DNAPairReadGraphBuilder::addInfoAbout2Read(std::string readName, int target, seqan::BamAlignmentRecord read) {
            PairReadGraphBuilder::addInfoAbout2Read(readName, target, read);
            TRACE("addInfoAboutRead2 readName=" << readName << " v=" << target);
            read2DistToEnd[readName] = readDist(read);
        }*/

        void DNAPairReadGraphBuilder::incEdgeWeight(seqan::BamAlignmentRecord read1, seqan::BamAlignmentRecord read2) {
            TRACE("incEdgeWeight");
            std::string readName = cutReadName(read1);
            if (read1DistToEnd[readName] + read2DistToEnd[readName] >
                distBetweenPairReads) {
                return;
            }
            PairReadGraphBuilder::incEdgeWeight(read1, read2);
        }

        std::string DNAPairReadGraphBuilder::getLibColor() {
            TRACE("getLibColor");
            int color[3] = {rand() % 100, 255, rand() % 100};
            return colorToString(color);
        }

        ContigGraph::Lib::Type DNAPairReadGraphBuilder::getLibType() {
            TRACE("getLibColor");
            return ContigGraph::Lib::DNA_PAIR;
        }
    }
}
