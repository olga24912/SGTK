#include <fstream>
#include <seqan/bam_io.h>
#include "ExonsInfo.h"

namespace findExon {
    void ExonsInfo::addInfo(std::string contigName, int contigLen, seqan::BamAlignmentRecord read) {
        if (exons == nullptr || exons->getContigName() != contigName) {
            if (exons != nullptr) {
                exons->finish();
                exons->writeExonBlock(out);
                delete(exons);
            }

            cover.resize(contigLen + 1, 0);
            misCover.resize(contigLen + 1, 0);

            exons = new ExonInfo(contigName, contigLen, cover, misCover);
        }

        exons->addReadInfo(read);
    }

    void ExonsInfo::printInfo() {
        if (exons != nullptr) {
            exons->finish();
            exons->writeExonBlock(out);
            delete(exons);
        }
    }

    ExonsInfo::ExonsInfo(std::string fileName) : out(fileName) {}
}
