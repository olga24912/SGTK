//
// Created by olga on 16.10.16.
//

#include "ReadsSplitter.h"
#include "FastaToolsIn.h"
#include "FastaToolsOut.h"

void ReadsSplitter::findAndSplitNotAlignmentReads(string rnaReadsFileName, string alignmentFileName,
                                                  string resFileName1, string resFileName2) {
    cerr << "start split reads" << endl;
    unordered_set<string> alignmentReads = findAlignmentReads(alignmentFileName);

    FastaToolsIn ftin;
    ftin.parse(rnaReadsFileName);

    FastaToolsOut ftout1;
    ftout1.putFileName(resFileName1);

    FastaToolsOut ftout2;
    ftout2.putFileName(resFileName2);

    cerr << "start rewrite reads" << endl;

    while (ftin.next()) {
        string readName = ftin.currentName();

        if (alignmentReads.count(readName) == 0) {
            string seq = ftin.currentRef();

            reads[readName] = seq;

            int len = (int) seq.size() / 2;
            string readName1 = readName;
            readName1 += "/1";
            string readName2 = readName;
            readName2 += "/2";

            ftout1.write(readName1, seq.substr(0, len));
            ftout2.write(readName2, seq.substr(len, seq.size() - len));
        }
    }

    ftin.close();
    ftout1.close();
    ftout2.close();
}

unordered_set<string> ReadsSplitter::findAlignmentReads(string fileName) {
    cerr << "start read sam file" << endl;
    unordered_set<string> usedReads;

    BamFileIn bamFile;
    open(bamFile, fileName.c_str());

    BamHeader samHeader;
    readHeader(samHeader, bamFile);

    BamAlignmentRecord read;
    while (!atEnd(bamFile)) {
        readRecord(read, bamFile);
        string name = string(toCString(read.qName));
        if (read.rID != -1) {
            usedReads.insert(name);
        }
    }
    close(bamFile);
    return usedReads;
}