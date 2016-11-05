//
// Created by olga on 05.11.16.
//

#include "ContigMerger.h"

void ContigMerger::evaluate(string contigsINFileName, string samReads1FileName, string samReads2FileName,
                            string contigOUTFileName, string samOutFileName, string contig1Name, string contig2Name) {
    this->contig1Name = contig1Name;
    this->contig2Name = contig2Name;

    contig1Val = findContig(contigsINFileName, contig1Name);
    contig2Val = findContig(contigsINFileName, contig2Name);
    cerr << "find contigs: " << contig1Val << " \n" << contig2Val << endl;

    mergeContigs();
    cerr << "\n newContig: " << contigVal << "\n";

    writeContig(contigOUTFileName);
    cerr << "out contig\n";
}

string ContigMerger::findContig(string fileIn, string name) {
    SeqFileIn seqFileIn(fileIn.c_str());
    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;

    readRecords(ids, seqs, seqFileIn);

    for (unsigned i = 0; i < length(ids); ++i) {
        string contigName = string(toCString(ids[i]));
        string seq = SeqanUtils::dna5ToString(toCString(seqs[i]), length(seqs[i]));

        if (contigName == name) {
            return seq;
        }
    }

    return "";
}

void ContigMerger::mergeContigs() {
    string ns = "";
    for (int i = 0; i < countN; ++i) {
        ns += 'N';
    }

    contigVal = contig1Val + ns + contig2Val;
}

void ContigMerger::writeContig(string fileName) {
    SeqFileOut out(fileName.c_str());

    StringSet<seqan::CharString> ids;
    appendValue(ids, CharString(contigName.c_str()));

    StringSet<seqan::CharString> seqs;
    appendValue(seqs, Dna5String(contigVal));

    writeRecords(out, ids, seqs);
}
