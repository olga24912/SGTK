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

    BamFileIn fileIn1(samReads1FileName.c_str());
    BamFileIn fileIn2(samReads2FileName.c_str());

    ofstream samOut;
    samOut.open(samOutFileName, std::ofstream::out);

    StringSet<CharString> contigNames;
    appendValue(contigNames, contigName);
    NameStoreCache<StringSet<CharString> > contigNamesCache(contigNames);
    BamIOContext<StringSet<CharString> > bamIOContext(contigNames, contigNamesCache);

    BamFileOut fileOut(context(fileIn1), samOut, Sam());
    fileOut.context = bamIOContext;

    writeHeader(fileOut);
    writeReads(fileOut, fileIn1, fileIn2);

    close(fileOut);
    close(fileIn1);
    close(fileIn2);
    samOut.close();
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

void ContigMerger::writeReads(BamFileOut &out, BamFileIn &in1, BamFileIn &in2) {
    BamAlignmentRecord record1;
    BamAlignmentRecord record2;
    while (!atEnd(in1)) {
        readRecord(record1, in1);
        readRecord(record2, in2);

        BamAlignmentRecord res1 = record1;
        res1.qName += "/1";
        res1.rID = 0;
        res1._buffer = "";

        BamAlignmentRecord res2 = record2;
        res2.qName += "/2";
        res2.rID = 0;
        res2.beginPos += contig1Val.size() + countN;


        writeRecord(out, res1);
        writeRecord(out, res2);
    }
}
