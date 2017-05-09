#include <cmath>
#include "DifStatistic.h"

void DifStatistic::calculateStatistic(Filter *filter, std::string coordFile, int libNum, int step, int mxVal) {
    InfoAboutContigsAlig aligInfo;
    aligInfo.parseCoordFile(filter, coordFile);

    int cnt_box = (int) (log(mxVal * 1.0) * 1.0/log(step) + 1);

    int cnt[cnt_box][7];
    for (int i = 0; i < cnt_box; ++i) {
        for (int j = 0; j < 7; ++j) {
            cnt[i][j] = 0;
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

        int bxnum = cnt_box - 1;
        if (smxW != 0) {
            double bigger = mxW * 1.0/smxW;
            bxnum = (int)(std::log(bigger)/std::log(step));
            if (bxnum >= cnt_box) {
                bxnum = cnt_box - 1;
            }
        }

        if (eval == -1) continue;
        ErrorType status = isCorrectEdge(aligInfo, filter, eval);
        cnt[bxnum][status]++;
    }


    for (int i = 0; i < cnt_box - 1; ++i) {
        printf("weight dif in: %d-%d: OK - %d; OVERLAP - %d; PART_ALIG - %d; BIG_DIST - %d;"
                       " WRONG_ORDER - %d; DIF_CHR - %d; NA - %d;\n",
               (int)pow(step, i), (int)pow(step, i + 1) - 1, cnt[i][0], cnt[i][1], cnt[i][2], cnt[i][3], cnt[i][4], cnt[i][5], cnt[i][6]);
    }
    printf("weight dif in: >%d: OK - %d; OVERLAP - %d; PART_ALIG - %d; BIG_DIST - %d;"
                   " WRONG_ORDER - %d; DIF_CHR - %d; NA - %d;\n",
           (int)pow(step, (cnt_box - 1)), cnt[cnt_box - 1][0], cnt[cnt_box - 1][1], cnt[cnt_box - 1][2], cnt[cnt_box - 1][3],
           cnt[cnt_box - 1][4], cnt[cnt_box - 1][5], cnt[cnt_box - 1][6]);
}
