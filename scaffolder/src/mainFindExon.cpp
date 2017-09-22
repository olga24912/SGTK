#include <Logger/log_writers.hpp>
#include <seqan/bam_io.h>
#include "FindExon/ExonInfo.h"
#include "FindExon/ExonsInfo.h"

//argv[1] = inFileName bamFormat
//argv[2] = outFileName


std::vector< std::pair<std::string, int> > readHeader(seqan::BamFileIn& bamFile) {
            INFO("start initHeader");
            typedef seqan::FormattedFileContext<seqan::BamFileIn, void>::Type TBamContext;

            seqan::BamHeader sam_header;
            readHeader(sam_header, bamFile);
            TBamContext const &bamContext = seqan::context(bamFile);
            size_t contig_num = seqan::length(contigNames(bamContext));

            DEBUG("Contig num=" << contig_num);

            std::vector<std::pair <std::string, int> > contigs;

            for (int i = 0; i < static_cast<int>(contig_num); ++i) {
                contigs.push_back(std::make_pair(
                        std::string(seqan::toCString(contigNames(bamContext)[i])),
                        contigLengths(bamContext)[i]));
            }

            INFO("finish initHeader");

            return contigs;
}

int main(int argc, char **argv) {
    using namespace findExon;
    logging::create_console_logger("/home/olga/bio-project/bio_scaffolder/scaffolder/src/log.properties");

    if (argc < 3) {
        ERROR("expect 2 args: <inFileName> <outFileName>, but find " << argc - 1);
    }

    std::string inFileName = argv[1];
    std::string outFileName = argv[2];

    ExonsInfo exons(outFileName);
    seqan::BamFileIn bamFile;

    if (!open(bamFile, inFileName.c_str())) {
        std::cerr << "could not open file";
        return 0;
    }

    std::vector<std::pair<std::string, int> > contigs = readHeader(bamFile);

    seqan::BamAlignmentRecord read;
    //int cnt = 10000000;

    while (!seqan::atEnd(bamFile) /*&& cnt > 0*/) {
        seqan::readRecord(read, bamFile);

        exons.addInfo(contigs[read.rID].first, contigs[read.rID].second, read); //read.beginPos, read.beginPos + seqan::getAlignmentLengthInRef(read)
        //--cnt;
    }

    exons.printInfo();
    return 0;
}
