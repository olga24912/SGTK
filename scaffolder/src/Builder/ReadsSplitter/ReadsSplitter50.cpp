#include "ReadsSplitter50.h"

void ReadsSplitter50::splitReads(std::string rnaUnmappedReadsFileName, std::string resFileName1, std::string resFileName2) {
    INFO("start splitReads rnaFileName=" << rnaUnmappedReadsFileName << " resFileName1=" << resFileName1 << " resFileName2=" << resFileName2);

    seqan::SeqFileIn seqFileIn(rnaUnmappedReadsFileName.c_str());
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;

    seqan::readRecords(ids, seqs, seqFileIn);

    seqan::SeqFileOut out1(resFileName1.c_str());
    seqan::SeqFileOut out2(resFileName2.c_str());

    for (unsigned i = 0; i < length(ids); ++i) {
        std::string readName = std::string(seqan::toCString(ids[i]));
        std::string seq = SeqanUtils::dna5ToString(seqan::toCString(seqs[i]), seqan::length(seqs[i]));
        
        int len = (int) seq.size() / 2;
        splitRead(readName, seq, len, out1, out2);
    }

    INFO("finish splitReads");
}

