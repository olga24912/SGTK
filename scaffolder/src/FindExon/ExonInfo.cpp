#include "ExonInfo.h"

namespace findExon {

    void findExon::ExonInfo::writeExonBlock(std::ostream &out) {
        out << contigName << " ";

        for (int i = 0; i < exons.size(); ++i) {
            out << exons[i].first << " " << exons[i].second << " ; ";
        }
        out << "\n";
    }

    void findExon::ExonInfo::addReadInfo(int b, int e) {
        if (exons.size() == 0) {
            exons.push_back(std::make_pair(b, e));
        } else if (exons[exons.size() - 1].second + 100 < b) {
            exons.push_back(std::make_pair(b, e));
        } else {
            exons[exons.size() - 1].first = std::min(exons[exons.size() - 1].first, b);
            exons[exons.size() - 1].second = std::max(exons[exons.size() - 1].second, e);
        }
    }
}