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

        int RuleExonBlocks::checkFirstCoord(ContigGraph *graph, int e) {
            int v = graph->getEdgeFrom(e);
            int cb = graph->getEdgeCoordB1(e);
            int ce = graph->getEdgeCoordE1(e);

            ContigGraph::Exon ex;
            ex.b = cb;
            ex.e = ce;

            std::vector<ContigGraph::Exon> exons = graph->getExons(v);
            int pos = std::upper_bound(exons.begin(), exons.end(), ex) - exons.begin();
            if (pos >= exons.size()) {
                return 0;
            }
            int lastId = exons[exons.size() - 1].id;
            if (exons[pos].id == lastId) {
                return 0;
            }

            if (exons[pos - 1].id == lastId) {
                return 0;
            }

            if (infoAlig.isCorrectEdge(graph, e) == statistics::InfoAboutContigsAlig::OK) {
                cnt[1] += 1;
            } else {
                cnt[0] += 1;
            }
            graph->delEdge(e);
            return 1;
        }

        void RuleExonBlocks::checkSecondCoord(ContigGraph *graph, int e) {
            int v = graph->getEdgeTo(e);
            int cb = graph->getEdgeCoordB2(e);
            int ce = graph->getEdgeCoordE2(e);

            ContigGraph::Exon ex;
            ex.b = cb;
            ex.e = ce;

            std::vector<ContigGraph::Exon> exons = graph->getExons(v);
            int pos = std::lower_bound(exons.begin(), exons.end(), ex) - exons.begin();
            if (pos == 0) { return; }
            int frId = exons[0].id;

            if (exons[pos - 1].id == frId) { return; }

            if (exons[pos].id == frId) { return; }

            if (infoAlig.isCorrectEdge(graph, e) == statistics::InfoAboutContigsAlig::OK) {
                cnt[1] += 1;

                std::cerr << v << " " << e << ": \n";
                std::cerr << "coordReads: " << cb << " " << ce << "\n";
                std::cerr << "exons: ";
                for (auto exn : exons) {
                    std::cerr << "( " << exn.b << " " << exn.e << " " << exn.id << " " << exn.cov << ") " ;
                }
                std::cerr << "\n";
            } else {
                cnt[0] += 1;
            }
            graph->delEdge(e);
        }
    }
}
