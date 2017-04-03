#include "GenTest.h"

void GenTest::genContig() {
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

std::string GenTest::rev_compl(std::string s) {
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

void GenTest::genTest(std::string refFileName, std::string outRefFileName) {
    this->refFileName = refFileName;
    this->outRefFileNema = outRefFileName;

    seqan::SeqFileIn seqFileIn(refFileName.c_str());
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;
    seqan::readRecords(ids, seqs, seqFileIn);
    dnaSeq = SeqanUtils::dna5ToString(seqan::toCString(seqs[0]), seqan::length(seqs[0]));

    genContig();
}