#include "RuleExonBlocks.h"
#include <algorithm>
#include <Filter/Statistics/StrandStatistic.h>

namespace filter {
    namespace scaffolder {
        using namespace contig_graph;
        void RuleExonBlocks::simplifyGraph(ContigGraph *graph) {
            //infoAlig.parseCoordFile(graph, "data/out.coords");
            //stat = new filter::statistics::StrandStatistic(graph, infoAlig, "data/ref.gff", "data/input.gff", "out");

            std::vector<int> vert = graph->getVertexList();
            for (int v : vert) {
                std::vector<int> edges = graph->getEdges(v);
                for (int e : edges) {
                    checkEdge(graph, e);
                }
            }
            //std::cerr << "del OK: " << cnt[1] << " del Wrong:" << cnt[0] << "\n";
        }

        void RuleExonBlocks::checkEdge(ContigGraph *graph, int e) {
            int res = checkFirstCoord(graph, e);
            if (res == 0) {
                checkSecondCoord(graph, e);
            }
        }

        int isGoodExon1(std::vector<ContigGraph::Exon> exons, ContigGraph::Exon ex) {
            int res = 0;
            if (exons.size() == 0) {
                return 0;
            }
            int pos = std::upper_bound(exons.begin(), exons.end(), ex) - exons.begin();
            if (pos < exons.size()) {
                if (exons[pos].b - 5 <= ex.b && ex.b <= exons[pos].e + 5) {
                    res =  (res | (1 << 1));
                }
            }
            if ((res & (1 << 1)) == 0 && pos > 0) {
                --pos;
                if (exons[pos].b - 5 <= ex.b && ex.b <= exons[pos].e + 5) {
                    res =  (res | (1 << 1));
                }
            }

            if ((res & (1 << 1)) == 0 && pos + 1 < exons.size()) {
                ++pos;
            }

            int lastId = exons[exons.size() - 1].id;

            if (exons[pos].id != lastId) {
                res |= 1;
            }

            return res;
        }

        int isGoodExon2(std::vector<ContigGraph::Exon> exons, ContigGraph::Exon ex) {
            int pos = std::lower_bound(exons.begin(), exons.end(), ex) - exons.begin();
            int res = 0;
            if (exons.size() == 0) {
                return 0;
            }
            if (pos < exons.size()) {
                if (exons[pos].b - 5 <= ex.b && ex.b <= exons[pos].e + 5) {
                    res =  (res | (1 << 1));
                }
            }
            if ((res & (1 << 1)) == 0 && pos > 0) {
                --pos;
                if (exons[pos].b - 5 <= ex.b && ex.b <= exons[pos].e + 5) {
                    res =  (res | (1 << 1));
                } else {
                    pos += 1;
                }
            }

            if ((res & (1 << 1)) == 0 && pos > 0) {
                --pos;
            }

            int frId = exons[0].id;
            if (exons[pos].id != frId) {
                res |= 1;
            }

            return res;
        }

        int cntConnection(ContigGraph *graph, int e) {
            int u = graph->getEdgeFrom(e);
            int v = graph->getEdgeTo(e);

            int cnt = 0;
            std::vector<int> edegs = graph->getEdges(u);
            for (int ee : edegs) {
                if (graph->getEdgeFrom(ee) == u && graph->getEdgeTo(ee) == v) {
                    ++cnt;
                }
            }
            return cnt;
        }

        int RuleExonBlocks::checkFirstCoord(ContigGraph *graph, int e) {
            int v = graph->getEdgeFrom(e);
            int cb = graph->getEdgeCoordB1(e);
            int ce = graph->getEdgeCoordE1(e);

            ContigGraph::Exon ex;
            ex.b = ce;
            ex.e = ce;

            int frres = isGoodExon1(graph->getExons(v, 1), ex);
            int scres = isGoodExon1(graph->getExons(v, 2), ex);

            if ((frres == 3 || scres == 3  || frres == 1 || scres == 1) && !(frres == 2 || scres == 2)) {
                //if (cntConnection(graph, e) == 1) {
               //     if (infoAlig.isCorrectEdge(graph, e) == alig_info::InfoAboutContigsAlig::OK) {
               //         cnt[1] += 1;
               //         stat->calculateStatisticForOneEdge(e);
               //     } else {
               //         cnt[0] += 1;
               //     }
               // }
                graph->delEdge(e);
                return 1;
            }
            return 0;
        }


        void RuleExonBlocks::checkSecondCoord(ContigGraph *graph, int e) {
            int v = graph->getEdgeTo(e);
            int cb = graph->getEdgeCoordB2(e);
            int ce = graph->getEdgeCoordE2(e);

            ContigGraph::Exon ex;
            ex.b = cb;
            ex.e = cb;

            int frres = isGoodExon2(graph->getExons(v, 1), ex);
            int scres = isGoodExon2(graph->getExons(v, 2), ex);

            if ((frres == 3 || scres == 3  || frres == 1 || scres == 1) && !(frres == 2 || scres == 2)) {
              //  if (cntConnection(graph, e) == 1) {
              //      if (infoAlig.isCorrectEdge(graph, e) == alig_info::InfoAboutContigsAlig::OK) {
              //          cnt[1] += 1;
              //      } else {
              //          cnt[0] += 1;
              //      }
              //  }
                graph->delEdge(e);
            }
        }
    }
}
