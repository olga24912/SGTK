#include "GeneAnnotationContigs.h"

namespace filter {
    namespace alig_info {
        GeneFinder::Gene GeneAnnotationContigs::findGene(ContigGraph::Vertex v, int coord) {
            if (v.name[v.name.size() - 1] != 'v') {
                GeneFinder::Gene gene = geneFinder.getGeneByCoord(v.name, coord);
                if (gene.b == -1) gene.strand = '!';

                return gene;
            } else {
                int nc = v.len - coord;
                std::string nn = v.name.substr(0, v.name.size() - 4);

                GeneFinder::Gene gene = geneFinder.getGeneByCoord(nn, nc);
                if (gene.b == -1) {
                    gene.strand = '!';
                    return gene;
                }

                gene.strand = ('+' + '-') - gene.strand;
                gene.b = v.len - gene.b;
                gene.e = v.len - gene.e;
                std::swap(gene.b, gene.e);

                for (int i = 0; i < (int)gene.exons.size(); ++i) {
                    gene.exons[i].strand = ('+' + '-') - gene.exons[i].strand;
                    gene.exons[i].b = v.len - gene.exons[i].b;
                    gene.exons[i].e = v.len - gene.exons[i].e;
                    std::swap(gene.exons[i].b, gene.exons[i].e);
                }

                return gene;
            }
        }
    }
}