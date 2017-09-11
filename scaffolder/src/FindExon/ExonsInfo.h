#ifndef SCAFFOLDER_EXONSINFO_H
#define SCAFFOLDER_EXONSINFO_H

#include<vector>
#include "ExonInfo.h"

namespace findExon {

    class ExonsInfo {
    private:
        std::vector<ExonInfo> exons;
    public:
        void addInfo(std::string exonName, int b, int e);

        void printInfo(std::string fileName);
    };
}

#endif //SCAFFOLDER_EXONSINFO_H
