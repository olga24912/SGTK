#include "CntRightConnectionStatistic.h"
#include "InfoAboutContigsAlig.h"

void CntRightConnectionStatistic::calculateStatistic(Filter *filter, std::string coordFile, int libNum) {
    InfoAboutContigsAlig aligInfo;
    aligInfo.parseCoordFile(filter, coordFile);

    int cnt[7] = {0, 0, 0, 0, 0, 0, 0};

    int n = filter->getVertexCount();
    for (int v = 0; v < n; ++v) {
        std::vector<int> edges = filter->getEdges(v);
        for (int e : edges) {
            if (filter->getEdgeLib(e) != libNum) continue;
            ErrorType status = isCorrectEdge(aligInfo, filter, e);
            cnt[status]++;
        }
    }


    printStatistic("", cnt);
}
