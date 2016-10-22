//
// Created by olga on 16.10.16.
//

#include "SplitReadsTest.h"

void SplitReadsTest::genTest(string fi, string fo, string rf, int readLen) {
    SeqFileIn seqFileIn(fi.c_str());
    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;
    readRecords(ids, seqs, seqFileIn);
    string fullref = SeqanUtils::dna5ToString(toCString(seqs[0]), length(seqs[0]));

    SeqFileOut outRead(rf.c_str());
    SeqFileOut outCon(fo.c_str());


    int spos = (int) (fullref.size() / 2);
    string ref = smallRef(fullref, 500, 1000, spos);
    for (int i = 0; i < ref.size() - readLen; ++i) {
        string readName = "read";
        readName += to_string(i);
        SeqanUtils::writeRec(outRead, readName, ref.substr(i, readLen));
    }

    SeqanUtils::writeRec(outCon, "contig1", fullref.substr(0, spos));
    SeqanUtils::writeRec(outCon, "contig2", fullref.substr(spos, fullref.size() - spos));
}

string SplitReadsTest::smallRef(string ref, int exonLen, int intronLen, int spos) {
    int len1 = intronLen/2, len2 = intronLen - intronLen/2;

    return ref.substr(spos - len1 - exonLen, exonLen).append(ref.substr(spos + len2, exonLen));
}

