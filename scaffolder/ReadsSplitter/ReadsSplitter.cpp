//
// Created by olga on 16.10.16.
//

#include "ReadsSplitter.h"


void ReadsSplitter::splitRead(string readName, string seq, int len, SeqFileOut& out1, SeqFileOut& out2) {
    string readName1 = readName;
    readName1 += "/1";
    string readName2 = readName;
    readName2 += "/2";

    StringSet<seqan::CharString> ids;
    appendValue(ids, CharString(readName1.c_str()));

    StringSet<seqan::CharString> seqs;
    appendValue(seqs, Dna5String(seq.substr(0, len)));

    writeRecords(out1, ids, seqs);
    clear(ids);
    clear(seqs);

    appendValue(ids, CharString(readName2.c_str()));
    appendValue(seqs, Dna5String(seq.substr(len, seq.size() - len)));

    writeRecords(out2, ids, seqs);
}
