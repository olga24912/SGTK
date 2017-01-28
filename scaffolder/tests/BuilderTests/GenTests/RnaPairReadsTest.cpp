#include "RnaPairReadsTest.h"

void RnaPairReadsTest::genTest(std::string refFileName, std::string read1FileName, std::string read2FileName,
                               std::string outRefFileName) {
    this->refFileName = refFileName;
    this->read1FileName = read1FileName;
    this->read2FileName = read2FileName;
    this->outRefFileNema = outRefFileName;

    seqan::SeqFileIn seqFileIn(refFileName.c_str());
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;
    seqan::readRecords(ids, seqs, seqFileIn);
    dnaSeq = SeqanUtils::dna5ToString(seqan::toCString(seqs[0]), seqan::length(seqs[0]));

    genRnaSeq();
    genPairRead();
    genContig();
}

void RnaPairReadsTest::genRnaSeq() {
    for (int itr = 0; itr < 5; ++itr) {
        for (int x = itr * 2000; x < itr * 2000 + 1500; ++x) {
            rnaSeq += dnaSeq[x];
        }
    }
}

void RnaPairReadsTest::genPairRead() {
    seqan::SeqFileOut outRead1(read1FileName.c_str());
    seqan::SeqFileOut outRead2(read2FileName.c_str());

    for (int i = 0; i < (int)rnaSeq.size() - 400; i += 10) {
        std::string name1;
        std::string name2;

        std::stringstream ss1, ss2;
        ss1 << "read" << i << "_1";
        name1 = std::string(ss1.str());
        ss2 << "read" << i << "_2";
        name2 = std::string(ss2.str());

        std::string read1;
        for (int j = i; j < i + 100; ++j) {
            read1 += rnaSeq[j];
        }

        std::string read2;

        for (int j = i + 300; j < i + 400; ++j) {
            read2 += rnaSeq[j];
        }

        read2 = rev_compl(read2);

        SeqanUtils::writeRec(outRead1, name1, read1);
        SeqanUtils::writeRec(outRead2, name2, read2);
    }
}

void RnaPairReadsTest::genContig() {
    seqan::SeqFileOut outCont(outRefFileNema.c_str());

    for (int i = 0; i < 10; ++i) {
        std::stringstream ss;
        ss << "contig" << i;
        if (i % 2 == 0) {
            SeqanUtils::writeRec(outCont, std::string(ss.str()), dnaSeq.substr((unsigned) i * 1000, 1000));
        } else {
            SeqanUtils::writeRec(outCont, std::string(ss.str()), rev_compl(dnaSeq.substr((unsigned) i * 1000, 1000)));
        }
    }
}

std::string RnaPairReadsTest::rev_compl(std::string s) {
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

