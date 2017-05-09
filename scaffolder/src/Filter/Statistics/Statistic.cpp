#include <iostream>
#include "Statistic.h"

Statistic::ErrorType Statistic::isCorrectEdge(InfoAboutContigsAlig &aligInfo, Filter *filter, int e) {
    int v = filter->getEdgeFrom(e);
    int u = filter->getEdgeTo(e);
    std::vector<InfoAboutContigsAlig::Alignment> valig = aligInfo.getAlignment(v);
    std::vector<InfoAboutContigsAlig::Alignment> ualig = aligInfo.getAlignment(u);

    ErrorType status = ErrorType::NA;

    for (auto val : valig) {
        for (auto ual : ualig) {
            auto al1 = val;
            auto al2 = ual;
            if (al1.chrName[al1.chrName.size() - 1] == 'v') {
                std::swap(al1, al2);
            }

            if (al1.chrName != al2.chrName) {
                if (ErrorType::DIF_CHR < status) {
                    status = ErrorType::DIF_CHR;
                }
                continue;
            }

            if (al1.coordBegin > al2.coordBegin) {
                if (ErrorType::WRONG_ORDER < status) {
                    status = ErrorType::WRONG_ORDER;
                }
                continue;
            }

            if (al2.coordBegin - al1.coordBegin > MAX_DIST) {
                if (ErrorType::BIG_DIST < status) {
                    status = ErrorType::BIG_DIST;
                }
                continue;
            }

            if ((al1.coordEnd - al1.coordBegin) * 10 < 7 * filter->getTargetLen(v) &&
                (al2.coordEnd - al2.coordBegin) * 10 < 7 * filter->getTargetLen(u)) {
                if (ErrorType::PART_ALIG < status) {
                    status = ErrorType::PART_ALIG;
                }
                continue;
            }

            if ((al1.coordEnd - al2.coordBegin) > MIN_OVERLAP)  {
                if (ErrorType::OVERLAP < status) {
                    status = ErrorType::OVERLAP;
                }
                continue;
            }

            status = ErrorType::OK;
        }
    }

    return status;
}
