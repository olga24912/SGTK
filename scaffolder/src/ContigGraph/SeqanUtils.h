#ifndef SCAFFOLDER_SEQANUTILS_H
#define SCAFFOLDER_SEQANUTILS_H

#include <string>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>

//handling info from seqans classes
class SeqanUtils {
public:
    //return reads name without "/1", "/2" end.
    static std::string cutReadName(seqan::BamAlignmentRecord read);
    //translate dna5 to string
    static std::string dna5ToString(seqan::Dna5* seq, int len);
    //write record(name, seq) in fasta/fastq format to out
    static void writeRec(seqan::SeqFileOut& out, std::string name, std::string seq);
};


#endif //SCAFFOLDER_SEQANUTILS_H
