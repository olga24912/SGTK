#ifndef SCAFFOLDER_EXONINFO_H
#define SCAFFOLDER_EXONINFO_H

#include <string>
#include <vector>
#include <iostream>

namespace findExon {

    class ExonInfo {
    private:
        struct Exon {
            int b;
            int e;
            int genId;
            double cover;
        };

        std::string contigName;
        int len;
        std::vector<Exon> exons;

        std::vector<int>& cover;
        std::vector<int>& misCover;
    public:
        ExonInfo(std::string contigName, int len, std::vector<int>&cover, std::vector<int>&misCover) :
                contigName(contigName), len(len), cover(cover), misCover(misCover) {}

        std::string getContigName() {
            return contigName;
        };

        void addReadInfo(seqan::BamAlignmentRecord read);

        void writeExonBlock(std::ostream &out);

        void finish();

        void findExonCon();

        bool hasConnection(Exon &exon, Exon &exon1, std::vector<int> &next);
    };
}


#endif //SCAFFOLDER_EXONINFO_H
