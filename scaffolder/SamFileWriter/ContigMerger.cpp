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

    BamFileIn fileIn(samReads1FileName.c_str());

    ofstream samOut;
    samOut.open(samOutFileName, std::ofstream::out);
    BamFileOut fileOut(context(fileIn), samOut, Sam());

    writeHeader(fileOut);

    close(fileOut);
    samOut.close();

    close(fileIn);
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

void ContigMerger::writeHeader(BamFileOut& out) {
    BamHeader header;
    std::stringstream ss;

    resize(header, 2);

    header[0].type = seqan::BAM_HEADER_FIRST;
    resize(header[0].tags, 1);
    header[0].tags[0].i1 = "VN";
    header[0].tags[0].i2 = "1.4";

    header[1].type = seqan::BAM_HEADER_REFERENCE;
    resize(header[1].tags, 2);

    header[1].tags[0].i1 = "SN";
    header[1].tags[0].i2 = contigName;
    header[1].tags[1].i1 = "LN";
    ss << length(contigVal);
    header[1].tags[1].i2 = ss.str();

    seqan::writeHeader(out, header);
}
