#ifndef SCAFFOLDER_GENEANNOTATIONCONTIGS_H
#define SCAFFOLDER_GENEANNOTATIONCONTIGS_H

#include "GeneFinder.h"
#include "InfoAboutContigsAlig.h"

namespace filter {
    namespace alig_info {
        class GeneAnnotationContigs {
        private:
            GeneFinder geneFinder;
        public:
            GeneAnnotationContigs(GeneFinder gf): geneFinder (gf) {}
            GeneFinder::Gene findGene(ContigGraph::Vertex v, int coord);
        };


    }
}

#endif //SCAFFOLDER_GENEANNOTATIONCONTIGS_H
