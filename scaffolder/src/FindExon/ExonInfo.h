#ifndef SCAFFOLDER_EXONINFO_H
#define SCAFFOLDER_EXONINFO_H

#include <string>
#include <vector>

namespace findExon {

    class ExonInfo {
    private:
        std::string contigName;
        std::vector<std::pair<int, int> > exons;
    public:
        void addReadInfo();

        void writeExonBlock(std::ofstream &out);
    };
}


#endif //SCAFFOLDER_EXONINFO_H
