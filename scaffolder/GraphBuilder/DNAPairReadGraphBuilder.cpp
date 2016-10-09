//
// Created by olga on 08.10.16.
//

#include "DNAPairReadGraphBuilder.h"

void DNAPairReadGraphBuilder::setDistBetweenPairReads(int distBetweenPairReads) {
    DNAPairReadGraphBuilder::distBetweenPairReads = distBetweenPairReads;
}

void DNAPairReadGraphBuilder::addInfoAboutRead(string readName, int target, BamAlignmentRecord read) {
    PairReadGraphBuilder::addInfoAboutRead(readName, target, read);
    read1DistToEnd[readName] = readDist(read);
}

int DNAPairReadGraphBuilder::readDist(BamAlignmentRecord read) {
    if (!hasFlagRC(read)) {
        return (graph.getTargetLength(2 * read.rID) - read.beginPos);
    } else {
        return (read.beginPos + read.tLen);
    }
}

void DNAPairReadGraphBuilder::addInfoAbout2Read(string readName, int target, BamAlignmentRecord read) {
    PairReadGraphBuilder::addInfoAbout2Read(readName, target, read);
    read2DistToEnd[readName] = readDist(read);
}
