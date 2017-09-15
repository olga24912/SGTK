#ifndef SCAFFOLDER_EXONSINFO_H
#define SCAFFOLDER_EXONSINFO_H

#include<vector>
#include "ExonInfo.h"

namespace findExon {

    class ExonsInfo {
    private:
        std::vector<ExonInfo> exons;
        std::vector<int> cover;
        std::vector<int> misCover;
    public:
        void addInfo(std::string contigName, int contigLen, seqan::BamAlignmentRecord read);

        void printInfo(std::string fileName);
    };
}

#endif //SCAFFOLDER_EXONSINFO_H
