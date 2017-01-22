//
// Created by olga on 22.10.16.
//

#include "SplitterByUnmappedEnd.h"

void SplitterByUnmappedEnd::splitReads(string rnaFileName, string resFileName1, string resFileName2) {
    SeqFileOut out1(resFileName1.c_str());
    SeqFileOut out2(resFileName2.c_str());

    BamFileIn bamFile;
    open(bamFile, rnaFileName.c_str());

    BamHeader samHeader;
    readHeader(samHeader, bamFile);

    while (!atEnd(bamFile)) {
        BamAlignmentRecord read;
        readRecord(read, bamFile);

        string name =  string(toCString(read.qName));
        string seq = string(toCString(CharString(read.seq)));

        CigarElement<char, unsigned> * cigar = toCString(read.cigar);
        CigarElement<char, unsigned> first = cigar[0];
        CigarElement<char, unsigned> last = cigar[length(read.cigar) - 1];

        if (first.operation == 'S') {
            if (first.count > MIN_READ_LEN) {
                splitRead(name, seq, first.count, out1, out2);
            }
        } else if (last.operation == 'S') {
            if (last.count > MIN_READ_LEN) {
                splitRead(name, seq, seq.size() - last.count, out1, out2);
            }
        }
    }
}
