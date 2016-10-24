//
// Created by olga on 22.10.16.
//

#include "ReadsSplitter50.h"

void ReadsSplitter50::splitReads(string rnaUnmappedReadsFileName, string resFileName1, string resFileName2) {
    cerr << "start split reads" << endl;
    SeqFileIn seqFileIn(rnaUnmappedReadsFileName.c_str());
    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;

    readRecords(ids, seqs, seqFileIn);

    SeqFileOut out1(resFileName1.c_str());
    SeqFileOut out2(resFileName2.c_str());

    cerr << "start rewrite reads" << endl;

    for (unsigned i = 0; i < length(ids); ++i) {
        string readName = string(toCString(ids[i]));
        cerr << readName <<  " ";
        string seq = SeqanUtils::dna5ToString(toCString(seqs[i]), length(seqs[i]));
        cerr << seq <<" " << endl;
        reads[readName] = seq;

        int len = (int) seq.size() / 2;
        splitRead(readName, seq, len, out1, out2);
    }
}

