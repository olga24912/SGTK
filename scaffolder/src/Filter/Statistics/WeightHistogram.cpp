#include "WeightHistogram.h"

namespace filter {
    namespace statistics {
        void WeightHistogram::histogram(ContigGraph *graph, int lib, int step) {
            int mx = 0;
            int n = graph->getVertexCount();
            for (int v = 0; v < n; ++v) {
                std::vector<int> edges = graph->getEdges(v);
                for (int e : edges) {
                    if (graph->getEdgeLib(e) != lib) continue;
                    if (graph->getEdgeWeight(e) > mx) {
                        mx = graph->getEdgeWeight(e);
                    }
                }
            }

            printf("Max val = %d\n", mx);

            mx = 20;
            std::vector<int> cnt((int) mx / step + 1, 0);
            int cnt_edge = 0;
            int sum_wg = 0;
            for (int v = 0; v < n; ++v) {
                std::vector<int> edges = graph->getEdges(v);
                for (int e : edges) {
                    if (graph->getEdgeLib(e) != lib) continue;
                    if (graph->getEdgeWeight(e) < mx) cnt[graph->getEdgeWeight(e) / step]++;
                    ++cnt_edge;
                    sum_wg += graph->getEdgeWeight(e);
                }
            }

            printf("Avg val = %d\n", sum_wg / cnt_edge);
            printf("Cnt edge = %d\n", cnt_edge);
            //freopen("hist5.dat", "w", stdout);
            for (int i = 1; i <= mx / step; ++i) {
                printf("%d \t %d \t %f\%\n", i * step, cnt[i], cnt[i] * 100.0 / cnt_edge);
            }

            //system("gnuplot <<< \"set term png; plot \'hist.dat\' using 1:2 with line\" > hist.png");
            //system("gnome-open hist.png");
        }
    }
}
