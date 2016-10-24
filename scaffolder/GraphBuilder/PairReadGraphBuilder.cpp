//
// Created by olga on 08.10.16.
//

#include "PairReadGraphBuilder.h"
#include "../Tools/SeqanUtils.h"

int PairReadGraphBuilder::pairTarget(int id) {
    return id ^ 1;
}

void PairReadGraphBuilder::setFileName2(const string &fileName2) {
    PairReadGraphBuilder::fileName2 = fileName2;
}

void PairReadGraphBuilder::setFileName1(const string &fileName1) {
    PairReadGraphBuilder::fileName1 = fileName1;
}

void PairReadGraphBuilder::evaluate() {
    read1Target.clear();
    cerr << "START" << endl;
    handleReads();
    filterEdge();
    graph->setLibNum(libName);
}

void PairReadGraphBuilder::readHeaderInit() {
    typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;

    BamHeader sam_header;
    readHeader(sam_header, bamFile1);
    TBamContext const &bamContext = context(bamFile1);
    size_t contig_num = length(contigNames(bamContext));

    string name;
    for (int i = 0; i < static_cast<int>(contig_num); ++i) {
        int length = contigLengths(bamContext)[i];
        name = string(toCString(contigNames(bamContext)[i]));
        graph->addVertex(2*i, name, 0,  length);
        name += "-rev";
        graph->addVertex(2*i + 1, name, 0, length);
    }
}

pair<string,int> PairReadGraphBuilder::processOneFirstRead(BamAlignmentRecord read) {
    readRecord(read, bamFile1);
    string readName = SeqanUtils::cutReadName(read);

    assert(read1Target.count(readName) == 0);

    bool isRev = hasFlagRC(read);
    int target = 2 * (read.rID);
    if (target < 0) {
        return make_pair(readName, -1);
    }
    if (isRev) {
        target++;
    }

    addInfoAboutRead(readName, target, read);
    return make_pair(readName, target);
}

void PairReadGraphBuilder::addInfoAboutRead(string readName, int target, BamAlignmentRecord read) {
    read1Target[readName] = target;
    addInfoAboutCover(target, read);
}

void PairReadGraphBuilder::addInfoAbout2Read(string readName, int target, BamAlignmentRecord read) {
    read2Target[readName] = target;
    addInfoAboutCover(target, read);
}

void PairReadGraphBuilder::addInfoAboutCover(int target, const BamAlignmentRecord &read) {
    int readLength = getAlignmentLengthInRef(read);
    int contigLength = graph->getTargetLength(target);
    graph->incTargetCover(target, static_cast<double>(readLength) / contigLength);
    graph->incTargetCover(target, static_cast<double>(readLength) / contigLength);
}


pair<string, int> PairReadGraphBuilder::processOneSecondRead(BamAlignmentRecord read) {
    readRecord(read, bamFile2);
    string readName = SeqanUtils::cutReadName(read);

    bool isRev = hasFlagRC(read);
    int target = 2 * (read.rID);
    if (target < 0) {
        return make_pair("", -1);
    }
    if ((!isRev && !oneSideRead) || (isRev && oneSideRead)) {
        target++;
    }

    addInfoAbout2Read(readName, target, read);
    return make_pair(readName, target);
}

void PairReadGraphBuilder::filterEdge() {
    graph->filterByContigLen(minContigLen);
    graph->filterByEdgeWeight(minEdgeWight);
}

void PairReadGraphBuilder::incEdgeWeight(string readName, int target) {
    if (read1Target.count(readName)) {
        if (read1Target[readName] == target ||
            read1Target[readName] == pairTarget(target)) {
            return;
        }
        int verFID = read1Target[readName], verSID = target,
                verRFID = pairTarget(verFID), verRSID = pairTarget(verSID);
        graph->incEdgeWeight(verFID, verSID);
        graph->incEdgeWeight(verRSID, verRFID);
    }
}

void PairReadGraphBuilder::setOneSideReadFlag(bool flag) {
    oneSideRead = flag;
}

void PairReadGraphBuilder::handleReads() {
    open(bamFile1, fileName1.c_str());
    open(bamFile2, fileName2.c_str());

    BamHeader samHeader2;
    readHeader(samHeader2, bamFile2);

    if (graph->getLibNum() == 1) {
        readHeaderInit();
    } else {
        BamHeader samHeader1;
        readHeader(samHeader1, bamFile1);
    }

    BamAlignmentRecord read1, read2;

    while (!atEnd(bamFile1) || !atEnd(bamFile2)) {
        pair<string, int> readInfo1 = make_pair("", -1);
        pair<string, int> readInfo2 = make_pair("", -1);

        if (!atEnd(bamFile1)) {
            readInfo1 = processOneFirstRead(read1);
        }

        if (!atEnd(bamFile2)) {
            readInfo2 = processOneSecondRead(read2);
        }
        if (readInfo2.second != -1 && read1Target.count(readInfo2.first)) {
            incEdgeWeight(readInfo2.first, readInfo2.second);
            read2Target.erase(readInfo2.first);
        }
        read1Target.erase(readInfo2.first);

        if (readInfo1.second != -1 && read2Target.count(readInfo1.first)) {
            incEdgeWeight(readInfo1.first, read2Target[readInfo1.first]);
            read1Target.erase(readInfo1.first);
        }
        read2Target.erase(readInfo1.first);
    }

    close(bamFile1);
    close(bamFile2);
}

