#include <Filter/AligInfo/GeneAnnotationContigs.h>
#include "StrandStatistic.h"

void
filter::statistics::StrandStatistic::calculateStatistic() {
    //std::ofstream outs(out);

    std::vector<int> vl = graph->getVertexList();
    for (int v : vl) {
        std::vector<int> edges = graph->getEdges(v);

        for (int i = 0; i < (int)(edges.size()); ++i) {
            int e = edges[i];
            int was = 0;
            for (int j = 0; j < i; ++j) {
                if (abs(graph->getEdgeCoordE1(e) - graph->getEdgeCoordE1(edges[j])) < 50) {
                    was = 1;
                }
            }
            if (was == 1) continue;
            if (graph->getEdgeWeight(e) < 5) continue;
            if (aligInfo.isCorrectEdge(graph, e) == alig_info::InfoAboutContigsAlig::OK) {
                calculateStatisticForOneEdge(e);
            }
        }
    }

    //outs.close();
}

void filter::statistics::StrandStatistic::printInfo(filter::contig_graph::ContigGraph *graph,
                                const filter::alig_info::InfoAboutContigsAlig &aligInfo, std::ofstream &outs, int v,
                                int e, const filter::alig_info::GeneFinder::Gene &gener,
                                const filter::alig_info::GeneFinder::Gene &genec) {
    printAligInfo(aligInfo, graph, e, outs);



    std::cerr << "printsmt";
    outs << e << " wieght="<< graph->getEdgeWeight(e) << " " << graph->getEdgeCoordB1(e) << " " << graph->getEdgeCoordE1(e) << " " << graph->getEdgeCoordB2(e) << " " << graph->getEdgeCoordE2(e) << "\n";
    outs << "\n";
    if (v == graph->getEdgeFrom(e)) {
        outs << "first part\n";
    } else {
        outs << "second part\n";
    }
    outs << v << " len=" << graph->getTargetLen(v) << " name=" << graph->getTargetName(v) << "\n";

    outs << "ref: \n";
    printGene(gener, outs);
    outs << "contigs: \n";
    printGene(genec, outs);
    outs << "---------------------\n" << std::flush;
}

void filter::statistics::StrandStatistic::printGene(filter::alig_info::GeneFinder::Gene gene, std::ofstream &outs) {
    if (gene.strand != '!') {
        outs << gene.strand << " " << gene.b << " " << gene.e << "\n";
        for (int i = 0; i < (int)gene.exons.size(); ++i) {
            outs << gene.exons[i].strand << " "<< gene.exons[i].b << " " << gene.exons[i].e << "\n";
        }
    }
}

bool filter::statistics::StrandStatistic::printAligInfo(filter::alig_info::InfoAboutContigsAlig alig,
                                                        filter::contig_graph::ContigGraph *graph, int e,
                                                        std::ofstream &outs) {
    std::vector<filter::alig_info::InfoAboutContigsAlig::Alignment> valig = alig.getAlignment(graph->getEdgeFrom(e));
    std::vector<filter::alig_info::InfoAboutContigsAlig::Alignment> ualig = alig.getAlignment(graph->getEdgeTo(e));

    int v = graph->getEdgeFrom(e);
    int u = graph->getEdgeTo(e);

    for (auto val : valig) {
        if ((val.coordEnd - val.coordBegin) * 1.0 / graph->getTargetLen(v) < 0.9) continue;
        for (auto ual : ualig) {
            if ((ual.coordEnd - ual.coordBegin) * 1.0 / graph->getTargetLen(u) < 0.9) continue;
            auto al1 = val;
            auto al2 = ual;

            if (al1.chrName != al2.chrName) {
                continue;
            }

            if (al1.coordBegin > al2.coordBegin) {
                continue;
            }

            if (al2.coordBegin - al1.coordEnd > 10000) {
                continue;
            }

            if ((al1.coordEnd - al2.coordBegin) > 500) {
                continue;
            }
            outs << al1.chrName << " " << al1.coordBegin << " " << al1.coordEnd << " " << al2.coordBegin << " "
                 << al2.coordEnd << "\n";
        }
    }
    return 1;
}

void filter::statistics::StrandStatistic::calculateStatisticForOneEdge(int e) {
    calculateStatisticForOneEdge1(e);
    calculateStatisticForOneEdge2(e);
}

void filter::statistics::StrandStatistic::calculateStatisticForOneEdge1(int e) {
    int v = graph->getEdgeFrom(e);
    filter::alig_info::GeneFinder::Gene gener = geneans.findGene(graph->getVertex(v), graph->getEdgeCoordE1(e));
    filter::alig_info::GeneFinder::Gene genec = geneContig.findGene(graph->getVertex(v), graph->getEdgeCoordE1(e));

    if (gener.strand != '+' && gener.strand != '-') {
        filter::statistics::StrandStatistic::printInfo(graph, aligInfo, outnone, v, e, gener, genec);
    }

    if (genec.strand != '+' && genec.strand != '-') {
        filter::statistics::StrandStatistic::printInfo(graph, aligInfo, outnone, v, e, gener, genec);
    } else if (gener.strand == genec.strand) {
        int cord = graph->getEdgeCoordE1(e);

        if (genec.b - 10 <= cord && cord <= genec.e + 10 ) {
            filter::statistics::StrandStatistic::printInfo(graph, aligInfo, outrht, v, e, gener, genec);
        } else {
            filter::statistics::StrandStatistic::printInfo(graph, aligInfo, outrhtwng, v, e, gener, genec);
        }
    } else {
        int cord = graph->getEdgeCoordE1(e);

        if (genec.b - 10 <= cord && cord <= genec.e + 10) {
            filter::statistics::StrandStatistic::printInfo(graph, aligInfo, outwngrht, v, e, gener, genec);
        } else {
            filter::statistics::StrandStatistic::printInfo(graph, aligInfo, outwng, v, e, gener, genec);
        }
    }
}

void filter::statistics::StrandStatistic::calculateStatisticForOneEdge2(int e) {
    int v = graph->getEdgeTo(e);
    filter::alig_info::GeneFinder::Gene gener = geneans.findGene(graph->getVertex(v), graph->getEdgeCoordB2(e));
    filter::alig_info::GeneFinder::Gene genec = geneContig.findGene(graph->getVertex(v), graph->getEdgeCoordB2(e));

    if (gener.strand != '+' && gener.strand != '-') {
        filter::statistics::StrandStatistic::printInfo(graph, aligInfo, outnone, v, e, gener, genec);
    }

    if (genec.strand != '+' && genec.strand != '-') {
        filter::statistics::StrandStatistic::printInfo(graph, aligInfo, outnone, v, e, gener, genec);
    } else if (gener.strand == genec.strand) {
        int cord = graph->getEdgeCoordB2(e);

        if (genec.b - 10 <= cord && cord <= genec.e + 10 ) {
            filter::statistics::StrandStatistic::printInfo(graph, aligInfo, outrht, v, e, gener, genec);
        } else {
            filter::statistics::StrandStatistic::printInfo(graph, aligInfo, outrhtwng, v, e, gener, genec);
        }
    } else {
        int cord = graph->getEdgeCoordB2(e);

        if (genec.b - 10 <= cord && cord <= genec.e + 10) {
            filter::statistics::StrandStatistic::printInfo(graph, aligInfo, outwngrht, v, e, gener, genec);
        } else {
            filter::statistics::StrandStatistic::printInfo(graph, aligInfo, outwng, v, e, gener, genec);
        }
    }
}
