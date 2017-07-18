#ifndef SCAFFOLDER_INFOABOUTCONTIGSALIG_H
#define SCAFFOLDER_INFOABOUTCONTIGSALIG_H


#include <Filter/Filters/ContigGraph.h>

namespace filter {
    namespace statistics {
        class InfoAboutContigsAlig {
        public:
            struct Alignment {
                std::string chrName;
                int coordBegin;
                int coordEnd;

                Alignment() {}

                Alignment(std::string chrName, int cb, int ce) : chrName(chrName), coordBegin(cb), coordEnd(ce) {}
            };

            void parseCoordFile(ContigGraph *graph, std::string coordFileName);

            std::vector<Alignment> getAlignment(int vertId);

            enum ErrorType {
                OK, OVERLAP, PART_ALIG, BIG_DIST, WRONG_ORDER, DIF_CHR, NA
            };
            const int MAX_DIST = 1000000;
            const int MIN_OVERLAP = 100;

            ErrorType isCorrectEdge(ContigGraph *filter, int e);

        private:
            std::vector<std::vector<Alignment> > alignment;
        };
    }
}

#endif //SCAFFOLDER_INFOABOUTCONTIGSALIG_H
