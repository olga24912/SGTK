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
    firstReads();
    cerr << "After first reads" << endl;
    secondReads();
    cerr << "After second reads" << endl;
    filterEdge();
}

void PairReadGraphBuilder::firstReads() {
    open(bamFile, fileName1.c_str());
    cerr << "first read" << endl;
    if (graph->getLibNum() == 0) {
        readHeaderInit();
    } else {
        BamHeader sam_header;
        readHeader(sam_header, bamFile);
    }
    BamAlignmentRecord read;
    while (!atEnd(bamFile)) {
        processOneFirstRead(read);
    }
    close(bamFile);
}

void PairReadGraphBuilder::readHeaderInit() {
    typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;

    BamHeader sam_header;
    readHeader(sam_header, bamFile);
    TBamContext const &bamContext = context(bamFile);
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

void PairReadGraphBuilder::processOneFirstRead(BamAlignmentRecord read) {
    readRecord(read, bamFile);
    string readName = SeqanUtils::cutReadName(read);

    cerr << "read: " << readName << "\n";

    assert(read1Target.count(readName) == 0);

    bool isRev = hasFlagRC(read);
    int target = 2 * (read.rID);
    if (target < 0) {
        return;
    }
    if (isRev) {
        target++;
    }

    addInfoAboutRead(readName, target, read);
}

void PairReadGraphBuilder::addInfoAboutRead(string readName, int target, BamAlignmentRecord read) {
    read1Target[readName] = target;
    addInfoAboutCover(target, read);
}

void PairReadGraphBuilder::addInfoAbout2Read(string readName, int target, BamAlignmentRecord read) {
    addInfoAboutCover(target, read);
}

void PairReadGraphBuilder::addInfoAboutCover(int target, const BamAlignmentRecord &read) {
    int readLength = getAlignmentLengthInRef(read);
    int contigLength = graph->getTargetLength(target);
    graph->incTargetCover(target, static_cast<double>(readLength) / contigLength);
    graph->incTargetCover(target, static_cast<double>(readLength) / contigLength);
}

void PairReadGraphBuilder::secondReads() {
    open(bamFile, fileName2.c_str());
    BamHeader samHeader;
    readHeader(samHeader, bamFile);
    BamAlignmentRecord read;
    pair<string, int> readInfo;
    while (!atEnd(bamFile)) {
        readInfo = processOneSecondRead(read);
        cerr << readInfo.first << " " << readInfo.second << endl;
        if (readInfo.second == -1) {
            continue;
        }
        incEdgeWeight(readInfo.first, readInfo.second);
    }
    close(bamFile);
}

pair<string, int> PairReadGraphBuilder::processOneSecondRead(BamAlignmentRecord read) {
    readRecord(read, bamFile);
    string readName = SeqanUtils::cutReadName(read);

    cerr << "read2: " << readName << endl;

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
        cerr << target << " " << verFID << endl;
        graph->incEdgeWeight(verFID, verSID);
        graph->incEdgeWeight(verRSID, verRFID);
    }
}

void PairReadGraphBuilder::setOneSideReadFlag(bool flag) {
    oneSideRead = flag;
}
