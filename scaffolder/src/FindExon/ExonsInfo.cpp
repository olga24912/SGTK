#include <fstream>
#include "ExonsInfo.h"

namespace findExon {
    void ExonsInfo::addInfo(std::string exonName, int b, int e) {
        if (exons.size() == 0 || exons[exons.size() - 1].getContigName() != exonName) {
            exons.push_back(ExonInfo(exonName));
        }

        exons[exons.size() - 1].addReadInfo(b, e);
    }

    void ExonsInfo::printInfo(std::string fileName) {
        std::ofstream out(fileName);

        for (auto &exon : exons) {
            exon.writeExonBlock(out);
        }
    }
}
