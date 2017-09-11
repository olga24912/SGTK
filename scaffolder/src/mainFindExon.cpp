#include <Logger/log_writers.hpp>
#include <seqan/bam_io.h>
#include "FindExon/ExonInfo.h"
#include "FindExon/ExonsInfo.h"

//argv[1] = inFileName bamFormat
//argv[2] = outFileName

int main(int argc, char **argv) {
    using namespace findExon;
    logging::create_console_logger("/home/olga/bio-project/bio_scaffolder/scaffolder/src/log.properties");

    if (argc < 3) {
        ERROR("expect 2 args: <inFileName> <outFileName>, but find " << argc - 1);
    }

    std::string inFileName = argv[1];
    std::string outFileName = argv[2];

    ExonsInfo exons;
    seqan::BamFileIn bamFile;

    if (!open(bamFile, inFileName.c_str())) {
        std::cerr << "could not open file";
        return 0;
    }


    return 0;
}
