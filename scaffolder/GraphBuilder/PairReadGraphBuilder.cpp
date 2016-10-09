//
// Created by olga on 08.10.16.
//

#include "PairReadGraphBuilder.h"

void PairReadGraphBuilder::setFileName2(const string &fileName2) {
    PairReadGraphBuilder::fileName2 = fileName2;
}

void PairReadGraphBuilder::setFileName1(const string &fileName1) {
    PairReadGraphBuilder::fileName1 = fileName1;
}

void PairReadGraphBuilder::evaluate() {
    cerr << "START" << endl;
    firstReads();
    cerr << "After first reads" << endl;
    secondReads();
    cerr << "After second reads" << endl;
}

void PairReadGraphBuilder::firstReads() {
    open(bamFile, fileName1.c_str());
    if (graph.getLibNum() == 0) {
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
        graph.addVertex(i, name, 0,  length);
        name += "-rev";
        graph.addVertex(i + 1, name, 0, length);
    }
}

void PairReadGraphBuilder::processOneFirstRead(BamAlignmentRecord read) {
    readRecord(read, bamFile);
    string readName = string(toCString(read.qName));

    if (readName.size() > 1) {
        if (readName[readName.size() - 2] == '/' &&
            readName[readName.size() - 1] == '1') {
            readName.resize(readName.size() - 2);
        }
    }

    assert(read1_target.count(read_name) == 0);

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
    auto readLength = getAlignmentLengthInRef(read);
    auto contigLength = getContigLength(read, bamFile);
    graph.incVertexCover(target, static_cast<double>(readLength) / contigLength);
    graph.incVertexCover(target, static_cast<double>(readLength) / contigLength);
}
