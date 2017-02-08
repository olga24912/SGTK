#include "ContigMerger.h"

void ContigMerger::evaluate(std::string contigsINFileName, std::string samReads1FileName, std::string samReads2FileName,
                            std::string contigOUTFileName, std::string samOutFileName,
                            std::string contig1Name, std::string contig2Name) {
    contig1Val = findContig(contigsINFileName, contig1Name);
    contig2Val = findContig(contigsINFileName, contig2Name);
    std::cerr << "find contigs: " << contig1Val << " \n" << contig2Val << std::endl;

    mergeContigs();
    std::cerr << "\n newContig: " << contigVal << "\n";

    writeContig(contigOUTFileName);
    std::cerr << "out contig\n";

    seqan::BamFileIn fileIn1(samReads1FileName.c_str());
    seqan::BamFileIn fileIn2(samReads2FileName.c_str());

    std::ofstream samOut;
    samOut.open(samOutFileName, std::ofstream::out);

    seqan::StringSet<seqan::CharString> contigNames;
    seqan::appendValue(contigNames, contigName);
    seqan::NameStoreCache<seqan::StringSet<seqan::CharString> > contigNamesCache(contigNames);
    seqan::BamIOContext<seqan::StringSet<seqan::CharString> > bamIOContext(contigNames, contigNamesCache);

    seqan::BamFileOut fileOut(seqan::context(fileIn1), samOut, seqan::Sam());
    fileOut.context = bamIOContext;

    writeHeader(fileOut);
    writeReads(fileOut, fileIn1, fileIn2);

    close(fileOut);
    close(fileIn1);
    close(fileIn2);
    samOut.close();
}

std::string ContigMerger::findContig(std::string fileIn, std::string name) {
    seqan::SeqFileIn seqFileIn(fileIn.c_str());
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;

    seqan::readRecords(ids, seqs, seqFileIn);

    for (unsigned i = 0; i < seqan::length(ids); ++i) {
        std::string contigName = getContigName(std::string(seqan::toCString(ids[i])));
        std::string seq = SeqanUtils::dna5ToString(seqan::toCString(seqs[i]), seqan::length(seqs[i]));

        std::cerr << contigName << " " << name << std::endl;

        if (contigName == name) {
            return seq;
        }
    }

    return "";
}

void ContigMerger::mergeContigs() {
    std::string ns = "";
    for (int i = 0; i < countN; ++i) {
        ns += 'N';
    }

    contigVal = contig1Val + ns + contig2Val;
}

void ContigMerger::writeContig(std::string fileName) {
    seqan::SeqFileOut out(fileName.c_str());

    seqan::StringSet<seqan::CharString> ids;
    seqan::appendValue(ids, seqan::CharString(contigName.c_str()));

    seqan::StringSet<seqan::CharString> seqs;
    seqan::appendValue(seqs, seqan::Dna5String(contigVal));

    seqan::writeRecords(out, ids, seqs);
}

void ContigMerger::writeHeader(seqan::BamFileOut& out) {
    std::stringstream ss;

    resize(header, 2);

    header[0].type = seqan::BAM_HEADER_FIRST;
    seqan::resize(header[0].tags, 1);
    header[0].tags[0].i1 = "VN";
    header[0].tags[0].i2 = "1.4";

    header[1].type = seqan::BAM_HEADER_REFERENCE;
    seqan::resize(header[1].tags, 2);

    header[1].tags[0].i1 = "SN";
    header[1].tags[0].i2 = contigName;
    header[1].tags[1].i1 = "LN";
    ss << seqan::length(contigVal);
    header[1].tags[1].i2 = ss.str();

    seqan::writeHeader(out, header);
}

void ContigMerger::writeReads(seqan::BamFileOut &out, seqan::BamFileIn &in1, seqan::BamFileIn &in2) {
    seqan::BamAlignmentRecord record1;
    seqan::BamAlignmentRecord record2;
    while (!seqan::atEnd(in1)) {
        seqan::readRecord(record1, in1);
        seqan::readRecord(record2, in2);

        seqan::BamAlignmentRecord res1 = record1;
        res1.qName += "/1";
        res1.rID = 0;
        res1._buffer = "";

        seqan::BamAlignmentRecord res2 = record2;
        res2.qName += "/2";
        res2.rID = 0;
        res2.beginPos += contig1Val.size() + countN;


        seqan::writeRecord(out, res1);
        seqan::writeRecord(out, res2);
    }
}

std::string ContigMerger::getContigName(std::string s) {
    std::string res = "";
    for (int i = 0; i < (int)s.size() && s[i] != ' '; ++i) {
        res += s[i];
    }
    return res;
}