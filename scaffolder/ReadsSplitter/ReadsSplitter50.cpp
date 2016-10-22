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

    FastaToolsOut ftout1;
    ftout1.putFileName(resFileName1);

    FastaToolsOut ftout2;
    ftout2.putFileName(resFileName2);

    cerr << "start rewrite reads" << endl;

    for (unsigned i = 0; i < length(ids); ++i) {
        string readName = string(toCString(ids[i]));
        string seq = SeqanUtils::dna5ToString(toCString(seqs[i]), length(seqs[i]));
        reads[readName] = seq;

        int len = (int) seq.size() / 2;
        splitRead(readName, seq, len, ftout1, ftout2);
    }

    ftout1.close();
    ftout2.close();
}

