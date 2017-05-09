#include "WeightStatistic.h"

void WeightStatistic::calculateStatistic(Filter *filter, std::string coordFile, int libNum, int step, int mxWeight) {
    InfoAboutContigsAlig aligInfo;
    aligInfo.parseCoordFile(filter, coordFile);

    int cnt_box = mxWeight/step + 1;

    int cnt[cnt_box][7];
    for (int i = 0; i < cnt_box; ++i) {
        for (int j = 0; j < 7; ++j) {
            cnt[i][j] = 0;
        }
    }

    int n = filter->getVertexCount();
    for (int v = 0; v < n; ++v) {
        std::vector<int> edges = filter->getEdges(v);
        for (int e : edges) {
            if (filter->getEdgeLib(e) != libNum) continue;
            ErrorType status = isCorrectEdge(aligInfo, filter, e);
            int wg = filter->getEdgeWeight(e)/step;
            if (wg >= cnt_box) {
                wg = cnt_box - 1;
            }
            cnt[wg][status]++;
        }
    }


    for (int i = 0; i < cnt_box - 1; ++i) {
        printf("Weight: %d-%d: OK - %d; OVERLAP - %d; PART_ALIG - %d; BIG_DIST - %d;"
                       " WRONG_ORDER - %d; DIF_CHR - %d; NA - %d;\n",
               i * step, (i + 1) * step - 1, cnt[i][0], cnt[i][1], cnt[i][2], cnt[i][3], cnt[i][4], cnt[i][5], cnt[i][6]);
    }
    printf("Weight: >%d: OK - %d; OVERLAP - %d; PART_ALIG - %d; BIG_DIST - %d;"
                   " WRONG_ORDER - %d; DIF_CHR - %d; NA - %d;\n",
           (cnt_box - 1)* step, cnt[cnt_box - 1][0], cnt[cnt_box - 1][1], cnt[cnt_box - 1][2], cnt[cnt_box - 1][3],
           cnt[cnt_box - 1][4], cnt[cnt_box - 1][5], cnt[cnt_box - 1][6]);
}
