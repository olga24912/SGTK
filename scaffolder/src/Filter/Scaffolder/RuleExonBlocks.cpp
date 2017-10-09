#include "RuleExonBlocks.h"
#include <algorithm>

namespace filter {
    namespace scaffolder {
        using namespace contig_graph;
        void RuleExonBlocks::simplifyGraph(ContigGraph *graph) {
            infoAlig.parseCoordFile(graph, "data/out.coords");

            std::vector<int> vert = graph->getVertexList();
            for (int v : vert) {
                std::vector<int> edges = graph->getEdges(v);
                for (int e : edges) {
                    checkEdge(graph, e);
                }
            }
            std::cerr << "del OK: " << cnt[1] << " del Wrong:" << cnt[0] << "\n";
        }

        void RuleExonBlocks::checkEdge(ContigGraph *graph, int e) {
            int res = checkFirstCoord(graph, e);
            if (res == 0) {
                checkSecondCoord(graph, e);
            }
        }

        bool isGoodExon1(std::vector<ContigGraph::Exon> exons, ContigGraph::Exon ex) {
            int pos = std::upper_bound(exons.begin(), exons.end(), ex) - exons.begin();
            if (pos >= exons.size()) {
                return true;
            }
            int lastId = exons[exons.size() - 1].id;

            if (exons[pos].id == lastId) {
                return true;
            }

            return false;
        }

        int RuleExonBlocks::checkFirstCoord(ContigGraph *graph, int e) {
            int v = graph->getEdgeFrom(e);
            int cb = graph->getEdgeCoordB1(e);
            int ce = graph->getEdgeCoordE1(e);

            ContigGraph::Exon ex;
            ex.b = cb;
            ex.e = ce;

            if (!isGoodExon1(graph->getExons(v, 1), ex) && !isGoodExon1(graph->getExons(v, 2), ex)) {
                if (infoAlig.isCorrectEdge(graph, e) == alig_info::InfoAboutContigsAlig::OK) {
                    cnt[1] += 1;
                } else {
                    cnt[0] += 1;
                }
                graph->delEdge(e);
                return 1;
            }
            return 0;
        }

        bool isGoodExon2(std::vector<ContigGraph::Exon> exons, ContigGraph::Exon ex) {
            int pos = std::lower_bound(exons.begin(), exons.end(), ex) - exons.begin();
            if (pos == 0) { return true; }
            int frId = exons[0].id;

            if (ex.b < exons[pos - 1].e && exons[pos - 1].id == frId) { return true; }

            if (exons[pos].id == frId) { return true; }

            return false;
        }

        void RuleExonBlocks::checkSecondCoord(ContigGraph *graph, int e) {
            int v = graph->getEdgeTo(e);
            int cb = graph->getEdgeCoordB2(e);
            int ce = graph->getEdgeCoordE2(e);

            ContigGraph::Exon ex;
            ex.b = cb;
            ex.e = ce;

            if (!isGoodExon2(graph->getExons(v, 1), ex) && !isGoodExon2(graph->getExons(v, 2), ex)) {

                if (infoAlig.isCorrectEdge(graph, e) == alig_info::InfoAboutContigsAlig::OK) {
                    cnt[1] += 1;
                } else {
                    cnt[0] += 1;
                }
                graph->delEdge(e);
            }
        }
    }
}
