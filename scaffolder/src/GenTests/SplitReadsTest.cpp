//
// Created by olga on 16.10.16.
//

#include "SplitReadsTest.h"

void SplitReadsTest::genTest(string fi, string fo, string rf, int readLen,int contigsCnt) {
    SeqFileIn seqFileIn(fi.c_str());
    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;
    readRecords(ids, seqs, seqFileIn);
    string fullref = SeqanUtils::dna5ToString(toCString(seqs[0]), length(seqs[0]));

    SeqFileOut outRead(rf.c_str());
    SeqFileOut outCon(fo.c_str());

    int delta = (int)(fullref.size()/contigsCnt);
    int p = 0;
    vector<int> spos(1, 0);
    for (int i = 0; i < contigsCnt - 1; ++i) {
        p += delta;
        spos.push_back(p);
    }
    spos.push_back((int)fullref.size());

    string ref = smallRef(fullref, 500, 1000, spos);
    for (int i = 0; i < ref.size() - readLen; ++i) {
        string readName = "read";
        readName += to_string(i);
        SeqanUtils::writeRec(outRead, readName, ref.substr(i, readLen));
    }

    for (int i = 0; i < contigsCnt; ++i) {
        string id;
        int x = i;
        while (x != 0) {
            id += '0' + x%10;
            x /= 10;
        }
        reverse(id.begin(), id.end());
        SeqanUtils::writeRec(outCon, "conitg" + id, fullref.substr(spos[i], spos[i + 1] - spos[i]));
    }
}

string SplitReadsTest::smallRef(string ref, int exonLen, int intronLen, vector<int> spos) {
    string res;
    int pos = intronLen;
    for (; pos + exonLen < ref.size(); pos += exonLen + intronLen) {
        res.append(ref.substr(pos, exonLen));
    }
    return res;
}

