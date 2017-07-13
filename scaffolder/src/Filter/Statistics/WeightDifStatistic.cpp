#include <cmath>
#include <iostream>
#include "WeightDifStatistic.h"

namespace filter {
    namespace statistics {
        void
        WeightDifStatistic::calculateStatistic(Filter *filter, std::string coordFile, int libNum, int step1, int mxW1,
                                               int step2, int mxW2) {
            InfoAboutContigsAlig aligInfo;
            aligInfo.parseCoordFile(filter, coordFile);

            int cnt_box1 = mxW1 / step1 + 1;
            int cnt_box2 = mxW2 / step2 + 1;

            int cnt[cnt_box1][cnt_box2][2][7];
            for (int i = 0; i < cnt_box1; ++i) {
                for (int j = 0; j < cnt_box2; ++j) {
                    for (int h = 0; h < 2; ++h) {
                        for (int g = 0; g < 7; ++g) {
                            cnt[i][j][h][g] = 0;
                        }
                    }
                }
            }

            int n = filter->getVertexCount();
            for (int v = 0; v < n; ++v) {
                std::vector<int> edges = filter->getEdges(v);
                int e1 = -1;
                int e2 = -1;
                for (int e : edges) {
                    if (filter->getEdgeLib(e) != libNum) {
                        continue;
                    }

                    if (e1 == -1 || filter->getEdgeWeight(e) > filter->getEdgeWeight(e1)) {
                        e2 = e1;
                        e1 = e;
                    } else if (e2 == -1 || filter->getEdgeWeight(e) > filter->getEdgeWeight(e2)) {
                        e2 = e;
                    }
                }

                if (e1 == -1) continue;
                InfoAboutContigsAlig::ErrorType status = aligInfo.isCorrectEdge(filter, e1);
                int wg1 = filter->getEdgeWeight(e1) / step1;
                int wg2 = 0;
                int flag = 0;
                if (e2 != -1) {
                    wg2 = filter->getEdgeWeight(e2) / step2;
                    if (wg2 >= cnt_box2) {
                        wg2 = cnt_box2 - 1;
                    }
                    if (filter->getEdgeTo(e1) == filter->getEdgeTo(e2)) {
                        flag = 1;
                    }
                }

                if (wg1 >= cnt_box1) {
                    wg1 = cnt_box1 - 1;
                }
                cnt[wg1][wg2][flag][status]++;
            }


            char str[100];
            for (int i = 0; i < cnt_box1 - 1; ++i) {
                for (int j = 0; j < cnt_box2 - 1; ++j) {
                    sprintf(str, "Weight 1: %d-%d \t weight 2: %d-%d;",
                            i * step1, (i + 1) * step1 - 1,
                            j * step2, (j + 1) * step2 - 1);
                    printStatistic(str, cnt[i][j][0]);
                }
                sprintf(str, "Weight 1: %d-%d \t weight 2: >%d;", i * step1, (i + 1) * step1 - 1,
                        (cnt_box2 - 1) * step2);
                printStatistic(str, cnt[i][cnt_box2 - 1][0]);
            }
            sprintf(str, "Weight 1: >%d \t weight 2: >%d;", (cnt_box1 - 1) * step1,
                    (cnt_box2 - 1) * step2);
            printStatistic(str, cnt[cnt_box1 - 1][cnt_box2 - 1][0]);
        }
    }
}
