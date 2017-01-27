#include "SplitReadsTest.h"

void SplitReadsTest::genTest(std::string fi, std::string fo, std::string rf, int readLen, int contigsCnt) {
    seqan::SeqFileIn seqFileIn(fi.c_str());
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;
    seqan::readRecords(ids, seqs, seqFileIn);
    std::string fullref = SeqanUtils::dna5ToString(seqan::toCString(seqs[0]), seqan::length(seqs[0]));

    seqan::SeqFileOut outRead(rf.c_str());
    seqan::SeqFileOut outCon(fo.c_str());

    int delta = (int)(fullref.size()/contigsCnt);
    int p = 0;
    std::vector<int> spos(1, 0);
    for (int i = 0; i < contigsCnt - 1; ++i) {
        p += delta;
        spos.push_back(p);
    }
    spos.push_back((int)fullref.size());

    std::string ref = smallRef(fullref, 500, 1000, spos);
    for (int i = 0; i < ref.size() - readLen; ++i) {
        std::string readName = "read";
        readName += std::to_string(i);
        SeqanUtils::writeRec(outRead, readName, ref.substr(i, readLen));
    }

    for (int i = 0; i < contigsCnt; ++i) {
        std::string id;
        int x = i;
        while (x != 0) {
            id += '0' + x%10;
            x /= 10;
        }
        std::reverse(id.begin(), id.end());
        SeqanUtils::writeRec(outCon, "conitg" + id, fullref.substr(spos[i], spos[i + 1] - spos[i]));
    }
}

std::string SplitReadsTest::smallRef(std::string ref, int exonLen, int intronLen, std::vector<int> spos) {
    std::string res;
    int pos = intronLen;
    for (; pos + exonLen < ref.size(); pos += exonLen + intronLen) {
        res.append(ref.substr(pos, exonLen));
    }
    return res;
}
