#include "TwoCompStatistic.h"

namespace filter {
    namespace statistics {
        void
        TwoCompStatistic::calculateStatistic(Filter *filter, std::string coordFile, int libNum1, int val1, int libNum2,
                                             int val2) {
            InfoAboutContigsAlig aligInfo;
            aligInfo.parseCoordFile(filter, coordFile);


            int cnt_box = 100;

            int cnt[cnt_box][7];
            for (int i = 0; i < cnt_box; ++i) {
                for (int g = 0; g < 7; ++g) {
                    cnt[i][g] = 0;
                }
            }


            int n = filter->getVertexCount();

            for (int v = 0; v < n; ++v) {
                std::vector<int> edges = filter->getEdges(v);
                int e1 = -1;
                int e1_2 = -1;
                int e2 = -1;
                for (int e : edges) {
                    if (filter->getEdgeLib(e) != libNum1) {
                        continue;
                    }

                    if (e1 == -1 || filter->getEdgeWeight(e) > filter->getEdgeWeight(e1)) {
                        e1 = e;
                    }
                }

                if (e1 == -1) continue;

                for (int e : edges) {
                    if (filter->getEdgeLib(e) != libNum1 || filter->getEdgeTo(e) == filter->getEdgeTo(e1)) {
                        continue;
                    }

                    if (e1_2 == -1 || filter->getEdgeWeight(e) > filter->getEdgeWeight(e1_2)) {
                        e1_2 = e;
                    }
                }


                for (int e : edges) {
                    if (filter->getEdgeLib(e) != libNum2 || filter->getEdgeTo(e) == filter->getEdgeTo(e1)) {
                        continue;
                    }

                    if (e2 == -1 || filter->getEdgeWeight(e) > filter->getEdgeWeight(e2)) {
                        e2 = e;
                    }
                }

                ErrorType status = isCorrectEdge(aligInfo, filter, e1);
                int wg1 = filter->getEdgeWeight(e1);
                if (wg1 >= cnt_box) wg1 = cnt_box - 1;
                int wg1_2 = 0;
                int wg2 = 0;

                if (e2 != -1) {
                    wg2 = filter->getEdgeWeight(e2);
                }

                if (e1_2 != -1) {
                    wg1_2 = filter->getEdgeWeight(e1_2);
                }

                if (wg1_2 != val1 || wg2 != val2) continue;
                cnt[wg1][status]++;
            }

            char str[100];

            for (int i = 0; i < cnt_box; ++i) {
                sprintf(str, "w: %d", i);
                printStatistic(str, cnt[i]);
            }
        }
    }
}