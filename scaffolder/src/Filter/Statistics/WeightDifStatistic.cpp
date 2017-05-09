#include <cmath>
#include <iostream>
#include "WeightDifStatistic.h"

void WeightDifStatistic::calculateStatistic(Filter *filter, std::string coordFile, int libNum, int stepW, int mxW,
                                            int stepDif, int mxDif) {
    InfoAboutContigsAlig aligInfo;
    aligInfo.parseCoordFile(filter, coordFile);

    int cnt_box_w = mxW/stepW + 1;
    int cnt_box_d = (int) (log(mxDif * 1.0) * 1.0/log(stepDif) + 1);

    int cnt[cnt_box_w][cnt_box_d][7];
    for (int i = 0; i < cnt_box_w; ++i) {
        for (int g = 0; g < cnt_box_d; ++g) {
            for (int j = 0; j < 7; ++j) {
                cnt[i][g][j] = 0;
            }
        }
    }

    int n = filter->getVertexCount();
    for (int v = 0; v < n; ++v) {
        std::vector<int> edges = filter->getEdges(v);
        int mxW = 0, smxW = 0;
        int eval = -1;
        for (int e : edges) {
            if (filter->getEdgeLib(e) != libNum) continue;
            if (filter->getEdgeWeight(e) > mxW) {
                smxW = mxW;
                mxW = filter->getEdgeWeight(e);
                eval = e;
            } else if (filter->getEdgeWeight(e) > smxW) {
                smxW = filter->getEdgeWeight(e);
            }
        }

        int bxDnum = cnt_box_d - 1;
        if (smxW != 0) {
            double bigger = mxW * 1.0/smxW;
            bxDnum = (int)(std::log(bigger)/std::log(stepDif));
            if (bxDnum >= cnt_box_d) {
                bxDnum = cnt_box_d - 1;
            }
        }

        if (eval == -1) continue;
        ErrorType status = isCorrectEdge(aligInfo, filter, eval);
        int bxWnum = filter->getEdgeWeight(eval)/stepW;
        if (bxWnum >= cnt_box_w) {
            bxWnum = cnt_box_w - 1;
        }
        cnt[bxWnum][bxDnum][status]++;
    }


    for (int g = 0; g < cnt_box_w - 1; ++g) {
        for (int i = 0; i < cnt_box_d - 1; ++i) {
            printf("weight: %d-%d; weight dif in: %d-%d: OK - %d; OVERLAP - %d; PART_ALIG - %d; BIG_DIST - %d;"
                           " WRONG_ORDER - %d; DIF_CHR - %d; NA - %d;\n",
                   g * stepW, (g + 1) * stepW - 1,
                   (int) pow(stepDif, i), (int) pow(stepDif, i + 1) - 1,
                   cnt[g][i][0], cnt[g][i][1], cnt[g][i][2], cnt[g][i][3],
                   cnt[g][i][4], cnt[g][i][5], cnt[g][i][6]);
        }
        printf("weight: %d-%d; weight dif in: >%d: OK - %d; OVERLAP - %d; PART_ALIG - %d; BIG_DIST - %d;"
                       " WRONG_ORDER - %d; DIF_CHR - %d; NA - %d;\n",
               g * stepW, (g + 1) * stepW - 1,
               (int) pow(stepDif, (cnt_box_d - 1)), cnt[g][cnt_box_d - 1][0],
               cnt[g][cnt_box_d - 1][1], cnt[g][cnt_box_d - 1][2],
               cnt[g][cnt_box_d - 1][3], cnt[g][cnt_box_d - 1][4],
               cnt[g][cnt_box_d - 1][5], cnt[g][cnt_box_d - 1][6]);
    }

    for (int i = 0; i < cnt_box_d - 1; ++i) {
        printf("weight: >%d; weight dif in: %d-%d: OK - %d; OVERLAP - %d; PART_ALIG - %d; BIG_DIST - %d;"
                       " WRONG_ORDER - %d; DIF_CHR - %d; NA - %d;\n",
               (cnt_box_w - 1) * stepW,
               (int) pow(stepDif, i), (int) pow(stepDif, i + 1) - 1,
               cnt[cnt_box_w - 1][i][0], cnt[cnt_box_w - 1][i][1], cnt[cnt_box_w - 1][i][2],
               cnt[cnt_box_w - 1][i][3], cnt[cnt_box_w - 1][i][4],
               cnt[cnt_box_w - 1][i][5], cnt[cnt_box_w - 1][i][6]);
    }
    printf("weight: >%d; weight dif in: >%d: OK - %d; OVERLAP - %d; PART_ALIG - %d; BIG_DIST - %d;"
                   " WRONG_ORDER - %d; DIF_CHR - %d; NA - %d;\n",
           (cnt_box_w - 1) * stepW,
           (int) pow(stepDif, (cnt_box_d - 1)),
           cnt[cnt_box_w - 1][cnt_box_d - 1][0],
           cnt[cnt_box_w - 1][cnt_box_d - 1][1], cnt[cnt_box_w - 1][cnt_box_d - 1][2],
           cnt[cnt_box_w - 1][cnt_box_d - 1][3], cnt[cnt_box_w - 1][cnt_box_d - 1][4],
           cnt[cnt_box_w - 1][cnt_box_d - 1][5], cnt[cnt_box_w - 1][cnt_box_d - 1][6]);
}
