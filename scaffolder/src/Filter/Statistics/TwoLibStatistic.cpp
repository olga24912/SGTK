#include "TwoLibStatistic.h"

namespace filter {
    namespace statistics {
        void
        TwoLibStatistic::calculateStatistic(Filter *filter, std::string coordFile, int libNum1, int step1, int mxW1,
                                            int step1_2, int mxW1_2, int libNum2, int step2, int mxW2) {
            InfoAboutContigsAlig aligInfo;
            aligInfo.parseCoordFile(filter, coordFile);

            int cnt_box1 = mxW1 / step1 + 1;
            int cnt_box1_2 = mxW1_2 / step1_2 + 1;
            int cnt_box2 = mxW2 / step2 + 1;

            int cnt[cnt_box1][cnt_box1_2][cnt_box2][7];
            for (int i = 0; i < cnt_box1; ++i) {
                for (int j = 0; j < cnt_box2; ++j) {
                    for (int h = 0; h < cnt_box1_2; ++h) {
                        for (int g = 0; g < 7; ++g) {
                            cnt[i][h][j][g] = 0;
                        }
                    }
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

                InfoAboutContigsAlig::ErrorType status = aligInfo.isCorrectEdge(filter, e1);
                int wg1 = filter->getEdgeWeight(e1) / step1;
                int wg1_2 = 0;
                int wg2 = 0;
                if (e2 != -1) {
                    wg2 = filter->getEdgeWeight(e2) / step2;
                    if (wg2 >= cnt_box2) {
                        wg2 = cnt_box2 - 1;
                    }
                }

                if (e1_2 != -1) {
                    wg1_2 = filter->getEdgeWeight(e1_2) / step1_2;
                    if (wg1_2 >= cnt_box1_2) {
                        wg1_2 = cnt_box1_2 - 1;
                    }
                }

                if (wg1 >= cnt_box1) {
                    wg1 = cnt_box1 - 1;
                }
                cnt[wg1][wg1_2][wg2][status]++;
            }

            char str[100];

            for (int i = 0; i < cnt_box1; ++i) {
                for (int g = 0; g < cnt_box1_2; ++g) {
                    for (int j = 0; j < cnt_box2; ++j) {
                        sprintf(str, "W_main: %d-%d \t w1: %d-%d \t w2: %d-%d;",
                                i * step1, (i + 1) * step1 - 1,
                                g * step1_2, (g + 1) * step1_2 - 1,
                                j * step2, (j + 1) * step2 - 1);
                        printStatistic(str, cnt[i][g][j]);
                    }
                }
            }
        }
    }
}
