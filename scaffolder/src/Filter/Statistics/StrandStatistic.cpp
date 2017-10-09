#include "StrandStatistic.h"

void
filter::statistics::StrandStatistic::calculateStatistic(filter::contig_graph::ContigGraph *graph, std::string coordFile,
                                                        std::string gffFile) {

    alig_info::InfoAboutContigsAlig aligInfo;
    aligInfo.parseCoordFile(graph, coordFile);
    alig_info::GeneAnnotationAns geneans(gffFile, aligInfo);

    std::vector<int> vl = graph->getVertexList();
    for (int v : vl) {
        std::vector<int> edges = graph->getEdges(v);

        for (int e : edges) {
            if (graph->getEdgeWeight(e) < 10) continue;
            if (aligInfo.isCorrectEdge(graph, e) == alig_info::InfoAboutContigsAlig::OK) {
                std::cerr << "-----\n";
                std::cerr << e << " "  << graph->getEdgeCoordB1(e) << " " << graph->getEdgeCoordE1(e) << " "
                          << graph->getEdgeCoordB2(e) << " " << graph->getEdgeCoordE2(e) << "\n";
                std::cerr << v << " len=" << graph->getTargetLen(v) << " name=" << graph->getTargetName(v) << "\n";
                alig_info::GeneFinder::Gene gene = geneans.findGene(graph->getVertex(v), graph->getEdgeCoordB1(e));
                if (gene.strand != '!') {
                    std::cerr << gene.strand << " " << gene.b << " " << gene.e << "\n";
                    for (int i = 0; i < (int)gene.exons.size(); ++i) {
                        std::cerr << gene.exons[i].strand << " "<< gene.exons[i].b << " " << gene.exons[i].e << "\n";
                    }
                }
                int u = graph->getEdgeTo(e);
                std::cerr << u << " len=" << graph->getTargetLen(u) << " name=" << graph->getTargetName(u) << "\n";
                gene = geneans.findGene(graph->getVertex(u), graph->getEdgeCoordB2(e));
                if (gene.strand != '!') {
                    std::cerr << gene.strand << " " << gene.b << " " << gene.e << "\n";
                    for (int i = 0; i < (int)gene.exons.size(); ++i) {
                        std::cerr << gene.exons[i].strand << " "<< gene.exons[i].b << " " << gene.exons[i].e << "\n";
                    }
                }
            }
        }

    }
}
