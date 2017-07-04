#include "RnaPairReadsTest.h"

void RnaPairReadsTest::genTest(std::string refFileName, std::string read1FileName, std::string read2FileName,
                               std::string outRefFileName) {
    this->refFileName = refFileName;
    this->read1FileName = read1FileName;
    this->read2FileName = read2FileName;
    this->outRefFileNema = outRefFileName;

    seqan::SeqFileIn seqFileIn(refFileName.c_str());
    //seqan::StringSet<seqan::CharString> ids;
    //seqan::StringSet<seqan::Dna5String> seqs;
    //seqan::readRecords(ids, seqs, seqFileIn);
    //dnaSeq = SeqanUtils::dna5ToString(seqan::toCString(seqs[0]), seqan::length(seqs[0]));

    genRnaSeq();
    genPairRead(rnaSeq);
    genContig();
}

void RnaPairReadsTest::genRnaSeq() {
    for (int itr = 0; itr < 5; ++itr) {
        for (int x = itr * 2000; x < itr * 2000 + 1500; ++x) {
            rnaSeq += dnaSeq[x];
        }
    }
}