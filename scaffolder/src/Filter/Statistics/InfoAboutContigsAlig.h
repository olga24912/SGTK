#ifndef SCAFFOLDER_INFOABOUTCONTIGSALIG_H
#define SCAFFOLDER_INFOABOUTCONTIGSALIG_H


#include <Filter/Filters/Filter.h>

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

            void parseCoordFile(Filter *graph, std::string coordFileName);

            std::vector<Alignment> getAlignment(int vertId);

        private:
            std::vector<std::vector<Alignment> > alignment;
        };
    }
}

#endif //SCAFFOLDER_INFOABOUTCONTIGSALIG_H
