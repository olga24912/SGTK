//
// Created by olga on 16.10.16.
//

#include "SeqanUtils.h"

string SeqanUtils::cutReadName(BamAlignmentRecord read) {
    string readName = string(toCString(read.qName));
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

string SeqanUtils::dna5ToString(Dna5 *seq, int len) {
    string res;

    string intToDNAChar = "ACGTN";

    for (int i = 0; i < len; ++i) {
        res += intToDNAChar[ordValue(seq[i])];
    }
    return res;
}

void SeqanUtils::writeRec(SeqFileOut& out, string name, string seq) {
    StringSet<seqan::CharString> ids;
    appendValue(ids, CharString(name.c_str()));

    StringSet<seqan::CharString> seqs;
    appendValue(seqs, Dna5String(seq));

    writeRecords(out, ids, seqs);
}
