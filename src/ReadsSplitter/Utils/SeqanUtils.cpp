#include "SeqanUtils.h"

namespace reads_splitter {
    namespace utils {
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