#include "SeqanUtils.h"

namespace builder {
    namespace contig_graph {
        std::string SeqanUtils::cutReadName(seqan::BamAlignmentRecord read) {
            std::string readName = std::string(seqan::toCString(read.qName));
            if (readName.size() > 1) {
                if ((readName[readName.size() - 2] == '/' ||
                     readName[readName.size() - 2] == '_') &&
                    (readName[readName.size() - 1] == '2' ||
                     readName[readName.size() - 1] == '1')) {
                    readName.resize(readName.size() - 2);
                }
            }
            return readName;
        }

        std::string SeqanUtils::dna5ToString(seqan::Dna5 *seq, int len) {
            std::string res;

            std::string intToDNAChar = "ACGTN";

            for (int i = 0; i < len; ++i) {
                res += intToDNAChar[ordValue(seq[i])];
            }
            return res;
        }

        void SeqanUtils::writeRec(seqan::SeqFileOut &out, std::string name, std::string seq) {
            seqan::StringSet<seqan::CharString> ids;
            seqan::appendValue(ids, seqan::CharString(name.c_str()));

            seqan::StringSet<seqan::CharString> seqs;
            seqan::appendValue(seqs, seqan::Dna5String(seq));

            seqan::writeRecords(out, ids, seqs);
        }
    }
}