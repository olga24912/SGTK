#include "GeneAnnotationAns.h"

namespace filter {
    namespace alig_info {
        GeneFinder::Gene GeneAnnotationAns::findGene(filter::contig_graph::ContigGraph::Vertex v, int coord) {
            INFO("start find gene for v=" << v.id << " " << coord);
            std::vector<InfoAboutContigsAlig::Alignment> aligns = aligInfo.getAlignment(v.id);

            for (InfoAboutContigsAlig::Alignment alig : aligns) {
                if ((alig.coordEnd - alig.coordBegin) * 1. / v.len > 0.9) {
                    if (alig.chrName[alig.chrName.size() - 1] == 'v') {
                        std::string nname = alig.chrName.substr(0, alig.chrName.size() - 4);
                        int nb = alig.chrLen - alig.coordEnd;
                        int ne = alig.chrLen - alig.coordBegin;
                        int nc = v.len - coord;

                        GeneFinder::Gene gene = geneFinder.getGeneByCoord(nname, nb + nc);
                        if (gene.b == -1) {
                            gene.strand = '!';
                            return gene;
                        }

                        gene.b = alig.chrLen - gene.b;
                        gene.e = alig.chrLen - gene.e;
                        std::swap(gene.b, gene.e);
                        gene.strand = ('+' + '-') - gene.strand;

                        for (int i = 0; i < gene.exons.size(); ++i) {
                            gene.exons[i].b = alig.chrLen - gene.exons[i].b;
                            gene.exons[i].e = alig.chrLen - gene.exons[i].e;
                            std::swap(gene.exons[i].b, gene.exons[i].e);
                            gene.exons[i].strand = ('+' + '-') - gene.exons[i].strand;
                        }

                        gene.b -= alig.coordBegin;
                        gene.e -= alig.coordBegin;

                        for (int i = 0; i < gene.exons.size(); ++i) {
                            gene.exons[i].b -= alig.coordBegin;
                            gene.exons[i].e -= alig.coordBegin;
                        }

                        return gene;
                    } else {
                        GeneFinder::Gene gene = geneFinder.getGeneByCoord(alig.chrName, alig.coordBegin + coord);
                        gene.b -= alig.coordBegin;
                        gene.e -= alig.coordBegin;

                        for (int i = 0; i < gene.exons.size(); ++i) {
                            gene.exons[i].b -= alig.coordBegin;
                            gene.exons[i].e -= alig.coordBegin;
                        }

                        return gene;
                    }
                }
            }

            GeneFinder::Gene gene;
            gene.strand = '!';
            return gene;
        }
    }
}
