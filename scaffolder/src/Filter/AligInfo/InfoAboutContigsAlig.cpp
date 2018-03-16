#include <map>
#include <fstream>
#include <ContigGraph/ContigGraph.h>
#include "InfoAboutContigsAlig.h"

namespace filter {
    namespace alig_info {
        using namespace contig_graph;
        void InfoAboutContigsAlig::parseCoordFile(ContigGraph *graph, std::string coordFileName) {
            if (alignment.size() > 0) return;
            int n = graph->getMaxVertId() + 1;
            alignment.resize(n);
            std::map<std::string, int> idByName;
            for (int i = 0; i < n; i += 2) {
                idByName[graph->getTargetName(i)] = i;
            }

            std::ifstream in(coordFileName);

            int l, r, lq, rq, x, chrLen;
            double xx;
            std::string rcont, qcont;

            while (in >> l >> r >> lq >> rq >> x >> x >> xx >> chrLen >> x >> rcont >> qcont) {
                int vid = idByName[qcont];
                if (lq > rq) {
                    vid ^= 1;
                    std::swap(lq, rq);
                }

                alignment[vid].push_back(Alignment(rcont, l, r, chrLen));
                alignment[vid ^ 1].push_back(Alignment(rcont + "-rev", chrLen - r, chrLen - l, chrLen));
            }

            in.close();
        }

        std::vector<InfoAboutContigsAlig::Alignment> InfoAboutContigsAlig::getAlignment(int vertId) {
            return alignment[vertId];
        }

        InfoAboutContigsAlig::ErrorType
        InfoAboutContigsAlig::isCorrectEdge(ContigGraph *filter, int e) {
            int v = filter->getEdgeFrom(e);
            int u = filter->getEdgeTo(e);
            std::vector<InfoAboutContigsAlig::Alignment> valig = getAlignment(v);
            std::vector<InfoAboutContigsAlig::Alignment> ualig = getAlignment(u);

            ErrorType status = ErrorType::NA;

            for (auto val : valig) {
                if ((val.coordEnd - val.coordBegin)*1.0/filter->getTargetLen(v) < 0.1) continue;
                for (auto ual : ualig) {
                    if ((ual.coordEnd - ual.coordBegin)*1.0/filter->getTargetLen(u) < 0.1) continue;
                    auto al1 = val;
                    auto al2 = ual;
                  /*  if (al1.chrName[al1.chrName.size() - 1] == 'v') {
                        std::swap(al1, al2);
                    }*/

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

                    if (al2.coordBegin - al1.coordEnd > MAX_DIST) {
                        if (ErrorType::BIG_DIST < status) {
                            status = ErrorType::BIG_DIST;
                        }
                        continue;
                    }

                    /*if ((al1.coordEnd - al1.coordBegin) * 10 < 7 * filter->getTargetLen(v) &&
                        (al2.coordEnd - al2.coordBegin) * 10 < 7 * filter->getTargetLen(u)) {
                        if (ErrorType::PART_ALIG < status) {
                            status = ErrorType::PART_ALIG;
                        }
                        continue;
                    }*/

                    if ((al1.coordEnd - al2.coordBegin) > MIN_OVERLAP) {
                        if (ErrorType::OVERLAP < status) {
                            status = ErrorType::OVERLAP;
                        }
                        continue;
                    }

                    status = ErrorType::OK;

                    DEBUG(al1.chrName << " " << al1.coordBegin << " " << al1.coordEnd << " " << al2.coordBegin << " " << al2.coordEnd << "\n");
                }
            }

            return status;
        }
    }
}