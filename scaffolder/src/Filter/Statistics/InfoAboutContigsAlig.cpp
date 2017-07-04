#include <map>
#include <fstream>
#include "InfoAboutContigsAlig.h"

namespace filter {
    namespace statistics {
        void InfoAboutContigsAlig::parseCoordFile(Filter *graph, std::string coordFileName) {
            int n = graph->getVertexCount();
            alignment.resize(n);
            std::map<std::string, int> idByName;
            for (int i = 0; i < n; i += 2) {
                idByName[graph->getTargetName(i)] = i;
            }

            std::ifstream in(coordFileName);

            int l, r, lq, rq, x;
            double xx;
            std::string rcont, qcont;

            while (in >> l >> r >> lq >> rq >> x >> x >> xx >> x >> x >> rcont >> qcont) {
                int vid = idByName[qcont];
                if (lq > rq) {
                    vid ^= 1;
                    std::swap(lq, rq);
                }

                alignment[vid].push_back(Alignment(rcont, l, r));
                alignment[vid ^ 1].push_back(Alignment(rcont + "-rev", l, r));
            }

            in.close();
        }

        std::vector<InfoAboutContigsAlig::Alignment> InfoAboutContigsAlig::getAlignment(int vertId) {
            return alignment[vertId];
        }
    }
}