#include "PairReadsTest.h"

void PairReadsTest::genPairRead(std::string ref) {
    seqan::SeqFileOut outRead1(read1FileName.c_str());
    seqan::SeqFileOut outRead2(read2FileName.c_str());

    for (int i = 0; i < (int)ref.size() - pairDist - 2*readLen; i += 10) {
        std::string name1;
        std::string name2;

        std::stringstream ss1, ss2;
        ss1 << "read" << i << "_1";
        name1 = std::string(ss1.str());
        ss2 << "read" << i << "_2";
        name2 = std::string(ss2.str());

        std::string read1;
        for (int j = i; j < i + readLen; ++j) {
            read1 += ref[j];
        }

        std::string read2;

        for (int j = i + readLen + pairDist; j < i + readLen * 2 + pairDist; ++j) {
            read2 += ref[j];
        }

        read2 = rev_compl(read2);

        SeqanUtils::writeRec(outRead1, name1, read1);
        SeqanUtils::writeRec(outRead2, name2, read2);
    }
}

void PairReadsTest::genContig() {
    seqan::SeqFileOut outCont(outRefFileNema.c_str());

    for (int i = 0; i < 10; ++i) {
        std::stringstream ss;
        ss << "contig" << i;
        if (i % 2 == 0) {
            SeqanUtils::writeRec(outCont, std::string(ss.str()), dnaSeq.substr((unsigned) i * contigLen, contigLen));
        } else {
            SeqanUtils::writeRec(outCont, std::string(ss.str()), rev_compl(dnaSeq.substr((unsigned) i * contigLen, contigLen)));
        }
    }
}

std::string PairReadsTest::rev_compl(std::string s) {
    std::string res = s;
    std::reverse(res.begin(), res.end());

    for (int i = 0; i < (int)res.size(); ++i) {
        if (res[i] == 'A') {
            res[i] = 'T';
        } else if (res[i] == 'T') {
            res[i] = 'A';
        } else if (res[i] == 'C') {
            res[i] = 'G';
        } else if (res[i] == 'G') {
            res[i] = 'C';
        }
    }

    return res;
}
