#ifndef SCAFFOLDER_EXONSINFO_H
#define SCAFFOLDER_EXONSINFO_H

#include<vector>
#include "ExonInfo.h"

namespace findExon {

    class ExonsInfo {
    private:
        ExonInfo* exons = nullptr;
        std::vector<int> cover;
        std::vector<int> misCover;

        std::ofstream out;
    public:
        ExonsInfo(std::string fileName);
        void addInfo(std::string contigName, int contigLen, seqan::BamAlignmentRecord read);

        void printInfo();
    };
}

#endif //SCAFFOLDER_EXONSINFO_H
