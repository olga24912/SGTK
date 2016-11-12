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

void PairReadGraphBuilder::setOneSideReadFlag(bool flag) {
    oneSideRead = flag;
}

void PairReadGraphBuilder::evaluate() {
    read1ByName.clear();
    cerr << "START" << endl;
    handleReads();
    filterEdge();
    graph->setLibName(libName);
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
    string readName = SeqanUtils::cutReadName(read);

    assert(read1ByName.count(readName) == 0);

    int target = get1Target(read);

    if (target < 0 || hasFlagSecondary(read)) {
        return make_pair(readName, -1);
    }
    addInfoAboutRead(readName, target, read);
    return make_pair(readName, target);
}

int PairReadGraphBuilder::get1Target(const BamAlignmentRecord &read) const {
    bool isRev = hasFlagRC(read);
    int target = 2 * (read.rID);
    if (isRev) {
        target++;
    }
    return target;
}

int PairReadGraphBuilder::get2Target(const BamAlignmentRecord &read) const {
    bool isRev = hasFlagRC(read);
    int target = 2 * (read.rID);
    if ((!isRev && !oneSideRead) || (isRev && oneSideRead)) {
        target++;
    }
    return target;
}

void PairReadGraphBuilder::addInfoAboutRead(string readName, int target, BamAlignmentRecord read) {
    read1ByName[readName] = read;
    addInfoAboutCover(target, read);
}

void PairReadGraphBuilder::addInfoAbout2Read(string readName, int target, BamAlignmentRecord read) {
    read2ByName[readName] = read;
    addInfoAboutCover(target, read);
}

void PairReadGraphBuilder::addInfoAboutCover(int target, const BamAlignmentRecord &read) {
    int readLength = getAlignmentLengthInRef(read);
    int contigLength = graph->getTargetLength(target);
    graph->incTargetCover(target, static_cast<double>(readLength) / contigLength);
    graph->incTargetCover(target, static_cast<double>(readLength) / contigLength);
}

pair<string, int> PairReadGraphBuilder::processOneSecondRead(BamAlignmentRecord read) {
    string readName = SeqanUtils::cutReadName(read);

    int target = get2Target(read);

    if (target < 0 || hasFlagSecondary(read)) {
        return make_pair("", -1);
    }
    addInfoAbout2Read(readName, target, read);
    return make_pair(readName, target);
}

void PairReadGraphBuilder::incEdgeWeight(BamAlignmentRecord read1, BamAlignmentRecord read2) {
    int target1 = get1Target(read1);
    int target2 = get2Target(read2);
    if (target1 == target2 || target1 == pairTarget(target2)) {
        return;
    }

    int verFID = target1, verSID = target2, verRFID = pairTarget(verFID), verRSID = pairTarget(verSID);
    int e1 = graph->incEdgeWeight(verFID, verSID);
    int e2 = graph->incEdgeWeight(verRSID, verRFID);

    samFileWriter.writeEdge(e1, read1, read2);
    samFileWriter.writeEdge(e2, read2, read1);
}

void PairReadGraphBuilder::handleReads() {
    open(bamFile1, fileName1.c_str());
    open(bamFile2, fileName2.c_str());

    samFileWriter.setFileIn(&bamFile1);

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
        pair<string, int> readInfo1;
        pair<string, int> readInfo2;

        if (!atEnd(bamFile1)) {
            readRecord(read1, bamFile1);
            readInfo1 = processOneFirstRead(read1);
        }
        if (!atEnd(bamFile2)) {
            readRecord(read2, bamFile2);
            readInfo2 = processOneSecondRead(read2);
        }
        if (readInfo2.first != "" && read1ByName.count(readInfo2.first)) {
            incEdgeWeight(read1ByName[readInfo2.first], read2);
            read2ByName.erase(readInfo2.first);
        }
        read1ByName.erase(readInfo2.first);

        if (readInfo1.first != "" && read2ByName.count(readInfo1.first)) {
            incEdgeWeight(read1, read2ByName[readInfo1.first]);
            read1ByName.erase(readInfo1.first);
        }
        read2ByName.erase(readInfo1.first);
    }

    close(bamFile1);
    close(bamFile2);
}

void PairReadGraphBuilder::filterEdge() {
    graph->filterByContigLen(minContigLen);
    graph->filterByEdgeWeight(minEdgeWight);
}