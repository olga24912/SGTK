#include "WeightStatistic.h"

namespace filter {
    namespace statistics {
        using namespace alig_info;
        void
        WeightStatistic::calculateStatistic(ContigGraph *graph, std::string coordFile, int libNum, int step, int mxWeight) {
            InfoAboutContigsAlig aligInfo;
            aligInfo.parseCoordFile(graph, coordFile);

            int cnt_box = mxWeight / step + 1;

            int cnt[cnt_box][7];
            for (int i = 0; i < cnt_box; ++i) {
                for (int j = 0; j < 7; ++j) {
                    cnt[i][j] = 0;
                }
            }

            int n = graph->getVertexCount();
            for (int v = 0; v < n; ++v) {
                std::vector<int> edges = graph->getEdges(v);
                for (int e : edges) {
                    if (graph->getEdgeLib(e) != libNum) continue;
                    InfoAboutContigsAlig::ErrorType status = aligInfo.isCorrectEdge(graph, e);
                    int wg = graph->getEdgeWeight(e) / step;
                    if (wg >= cnt_box) {
                        wg = cnt_box - 1;
                    }
                    cnt[wg][status]++;
                }
            }


            char str[100];
            for (int i = 0; i < cnt_box - 1; ++i) {
                sprintf(str, "Weight: %d-%d:",
                        i * step, (i + 1) * step - 1);
                printStatistic(str, cnt[i]);
            }

            sprintf(str, "Weight: >%d;", (cnt_box - 1) * step);
            printStatistic(str, cnt[cnt_box - 1]);
        }
    }
}
