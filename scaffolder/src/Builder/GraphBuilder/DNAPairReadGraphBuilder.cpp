#include "DNAPairReadGraphBuilder.h"
#include "Builder/Tools/SeqanUtils.h"

void DNAPairReadGraphBuilder::setDistBetweenPairReads(int distBetweenPairReads) {
    DNAPairReadGraphBuilder::distBetweenPairReads = distBetweenPairReads;
}

void DNAPairReadGraphBuilder::addInfoAboutRead(std::string readName, int target, seqan::BamAlignmentRecord read) {
    PairReadGraphBuilder::addInfoAboutRead(readName, target, read);
    read1DistToEnd[readName] = readDist(read);
}

int DNAPairReadGraphBuilder::readDist(seqan::BamAlignmentRecord read) {
    if (!seqan::hasFlagRC(read)) {
        return (graph->getTargetLength(2 * read.rID) - read.beginPos - (int)(read.seq.data_end - read.seq.data_begin));
    } else {
        return (read.beginPos);
    }
}

void DNAPairReadGraphBuilder::addInfoAbout2Read(std::string readName, int target, seqan::BamAlignmentRecord read) {
    PairReadGraphBuilder::addInfoAbout2Read(readName, target, read);
    read2DistToEnd[readName] = readDist(read);
}

void DNAPairReadGraphBuilder::incEdgeWeight(seqan::BamAlignmentRecord read1, seqan::BamAlignmentRecord read2) {
    std::string readName = SeqanUtils::cutReadName(read1);
    if (read1DistToEnd[readName] + read2DistToEnd[readName] >
        distBetweenPairReads) {
        return;
    }
    PairReadGraphBuilder::incEdgeWeight(read1, read2);
}

std::string DNAPairReadGraphBuilder::getLibColor() {
    int color[3] = {rand() % 100, 255, rand()%100};
    return colorToString(color);
}
