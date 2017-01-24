//
// Created by olga on 08.10.16.
//

#include "DNAPairReadGraphBuilder.h"
#include "Builder/Tools/SeqanUtils.h"

void DNAPairReadGraphBuilder::setDistBetweenPairReads(int distBetweenPairReads) {
    DNAPairReadGraphBuilder::distBetweenPairReads = distBetweenPairReads;
}

void DNAPairReadGraphBuilder::addInfoAboutRead(string readName, int target, BamAlignmentRecord read) {
    PairReadGraphBuilder::addInfoAboutRead(readName, target, read);
    read1DistToEnd[readName] = readDist(read);
}

int DNAPairReadGraphBuilder::readDist(BamAlignmentRecord read) {
    if (!hasFlagRC(read)) {
        return (graph->getTargetLength(2 * read.rID) - read.beginPos);
    } else {
        return (read.beginPos + read.tLen);
    }
}

void DNAPairReadGraphBuilder::addInfoAbout2Read(string readName, int target, BamAlignmentRecord read) {
    PairReadGraphBuilder::addInfoAbout2Read(readName, target, read);
    read2DistToEnd[readName] = readDist(read);
}

void DNAPairReadGraphBuilder::incEdgeWeight(BamAlignmentRecord read1, BamAlignmentRecord read2) {
    string readName = SeqanUtils::cutReadName(read1);
    if (read1DistToEnd[readName] + read2DistToEnd[readName] >
        distBetweenPairReads) {
        return;
    }
    PairReadGraphBuilder::incEdgeWeight(read1, read2);
}

string DNAPairReadGraphBuilder::getLibColor() {
    int color[3] = {rand() % 100, 255, rand()%100};
    return colorToString(color);
}
