#include <fstream>
#include <seqan/bam_io.h>
#include "ExonsInfo.h"

namespace findExon {
    void ExonsInfo::addInfo(std::string contigName, int contigLen, seqan::BamAlignmentRecord read) {
        if (exons.size() == 0 || exons[exons.size() - 1].getContigName() != contigName) {
            if (exons.size() > 0) {
                exons[exons.size() - 1].finish();
            }

            if (cover.size() <= contigLen) {
                cover.resize(contigLen + 1, 0);
                misCover.resize(contigLen + 1, 0);
            }

            exons.push_back(ExonInfo(contigName, contigLen, cover, misCover));
        }

        exons[exons.size() - 1].addReadInfo(read);
    }

    void ExonsInfo::printInfo(std::string fileName) {
        std::ofstream out(fileName);
        if (exons.size() > 0) {
            exons[exons.size() - 1].finish();
        }

        for (auto &exon : exons) {
            exon.writeExonBlock(out);
        }
    }
}
