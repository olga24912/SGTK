#include "SplitterByUnmappedEnd.h"

namespace reads_splitter {
    void
    SplitterByUnmappedEnd::splitReads(std::string rnaFileName, std::string resFileName1, std::string resFileName2) {
        INFO("start splitReads rnaFileName=" << rnaFileName << " resFileName1=" << resFileName1 << " resFileName2="
                                             << resFileName2);

        seqan::SeqFileOut out1(resFileName1.c_str());
        seqan::SeqFileOut out2(resFileName2.c_str());

        seqan::BamFileIn bamFile;
        seqan::open(bamFile, rnaFileName.c_str());

        seqan::BamHeader samHeader;
        seqan::readHeader(samHeader, bamFile);

        while (!seqan::atEnd(bamFile)) {
            seqan::BamAlignmentRecord read;
            seqan::readRecord(read, bamFile);

            std::string name = std::string(toCString(read.qName));
            std::string seq = std::string(toCString(seqan::CharString(read.seq)));

            seqan::CigarElement<char, unsigned> *cigar = toCString(read.cigar);
            seqan::CigarElement<char, unsigned> first = cigar[0];
            seqan::CigarElement<char, unsigned> last = cigar[length(read.cigar) - 1];

            if (first.operation == 'S') {
                if (first.count > MIN_READ_LEN) {
                    splitRead(name, seq, first.count, out1, out2);
                }
            } else if (last.operation == 'S') {
                if (last.count > MIN_READ_LEN) {
                    splitRead(name, seq, seq.size() - last.count, out1, out2);
                }
            }
        }

        INFO("finish splitReads");
    }
}
