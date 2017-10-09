#ifndef SCAFFOLDER_GENEANNOTATIONANS_H
#define SCAFFOLDER_GENEANNOTATIONANS_H


#include "GeneFinder.h"
#include "InfoAboutContigsAlig.h"

namespace filter {
    namespace alig_info {
        class GeneAnnotationAns {
        private:
            GeneFinder geneFinder;
            InfoAboutContigsAlig aligInfo;
        public:
            GeneAnnotationAns(GeneFinder gf, InfoAboutContigsAlig alig): geneFinder(gf), aligInfo(alig) {}

            GeneFinder::Gene findGene(ContigGraph::Vertex v, int coord);
        };
    }
}


#endif //SCAFFOLDER_GENEANNOTATIONANS_H
