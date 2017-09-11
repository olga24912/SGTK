#ifndef SCAFFOLDER_EXONINFO_H
#define SCAFFOLDER_EXONINFO_H

#include <string>
#include <vector>
#include <iostream>

namespace findExon {

    class ExonInfo {
    private:
        std::string contigName;
        std::vector<std::pair<int, int> > exons;
    public:
        ExonInfo(std::string contigName) : contigName(contigName) {}

        std::string getContigName() {
            return contigName;
        };

        void addReadInfo(int b1, int e1);

        void writeExonBlock(std::ostream &out);
    };
}


#endif //SCAFFOLDER_EXONINFO_H
