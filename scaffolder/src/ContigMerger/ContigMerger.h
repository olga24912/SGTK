#ifndef SCAFFOLDER_CONTIGMERGER_H
#define SCAFFOLDER_CONTIGMERGER_H

#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include "ContigGraph/SeqanUtils.h"
#include <string>
#include <Logger/logger.hpp>


namespace contig_merger {
//Merge two contig with some space between.
//Change alignment reads for this new contig
    class ContigMerger {
    private:
        seqan::BamHeader header;

        int countN = 300; //space between contigs in new contig

        std::string contig1Val;
        std::string contig2Val;

        std::string contigName = "contig"; //name of the result contig
        std::string contigVal; //value of the result contig

        std::string findContig(std::string fileIn, std::string name);

        void mergeContigs();

        void writeContig(std::string fileName);

        void writeHeader(seqan::BamFileOut &out);

        void writeReads(seqan::BamFileOut &out, seqan::BamFileIn &in1, seqan::BamFileIn &in2);

        std::string getContigName(std::string s);

    public:
        //merge two contigs and realignment reads
        //contgsInFileName - name of file with all contigs in fasta/fastq format
        //samReads1FileName - name of file with alignment reads on first contig in SAM/BAM format
        //samReads2FileName - name of file with alignment reads on second contig in SAM/BAM format
        //contigOUTFileName - name of file where will be written one merge contig in fasta/fsatq format
        //samOutFileName - name of file with alignment reads on new contig
        //contig1Name - name of  the first contig
        //contig2Name - name of the second contig
        void evaluate(std::string contigsINFileName, std::string samReads1FileName,
                      std::string samReads2FileName, std::string contigOUTFileName,
                      std::string samOutFileName,
                      std::string contig1Name, std::string contig2Name);

    private:
        DECL_LOGGER("ContigMerger");
    };
}

#endif //SCAFFOLDER_CONTIGMERGER_H
