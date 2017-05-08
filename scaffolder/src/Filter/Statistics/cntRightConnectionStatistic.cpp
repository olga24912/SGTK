#include "cntRightConnectionStatistic.h"
#include "InfoAboutContigsAlig.h"

void cntRightConnectionStatistic::calculateStatistic(Filter *filter, std::string coordFile, int libNum) {
    InfoAboutContigsAlig aligInfo;
    aligInfo.parseCoordFile(filter, coordFile);

    int cntRight = 0;
    int cntWrong = 0;
    int cntNA = 0;

    int n = filter->getVertexCount();
    for (int v = 0; v < n; ++v) {
        std::vector<int> edges = filter->getEdges(v);
        for (int e : edges) {
            if (filter->getEdgeLib(e) != libNum) continue;
            int u = filter->getEdgeTo(e);
            std::vector<InfoAboutContigsAlig::Alignment> valig = aligInfo.getAlignment(v);
            std::vector<InfoAboutContigsAlig::Alignment> ualig = aligInfo.getAlignment(u);

            int wasOk = 0;

            for (auto al1 : valig) {
                for (auto al2 : ualig) {
                    if (al1.coordEnd < al2.coordBegin && al1.chrName == al2.chrName &&
                            al2.coordEnd - al1.coordBegin < MAX_DIST) {
                        if ((al1.coordEnd - al1.coordBegin) * 10 > 9 * filter->getTargetLen(v) &&
                                (al2.coordEnd - al2.coordBegin) * 10 > 9 * filter->getTargetLen(u)) {
                            wasOk = 1;
                        }
                        if (wasOk != 1) wasOk = 2;
                    }
                }
            }

            if (wasOk == 1) {
                cntRight++;
            } else if (wasOk == 2){
                cntNA++;
            } else {
                cntWrong++;
            }
        }
    }

    printf("cntRightConnection: %d\n", cntRight);
    printf("cntWrongConnection: %d\n", cntWrong);
    printf("cntNAConnection: %d\n", cntNA);
}
