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
        name = contigNames(bamContext)[i];
        graph.addVertex(i, name, 0,  length);
        name += "-rev";
        graph.addVertex(i + 1, name, 0, length);
    }
}
