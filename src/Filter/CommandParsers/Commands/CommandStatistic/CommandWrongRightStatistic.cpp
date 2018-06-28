#include <set>
#include "CommandWrongRightStatistic.h"

namespace filter {
    namespace  commands {
        void filter::commands::CommandWrongRightStatistic::execute(std::string argv, filter::commands::State &state,
                                                                   filter::contig_graph::ContigGraph &graph) {
            std::stringstream ss(argv);
            std::string outFile;
            ss >> outFile;
            std::ofstream out(outFile);
            int wlib;
            std::vector<int> rlibs;

            ss >> wlib;
            int lb;
            while (ss >> lb) {
                rlibs.push_back(lb);
            }

            alig_info::InfoAboutContigsAlig aligInfo;
            aligInfo.parseCoordFile(&graph, state.coordFile);

            std::vector<int> vert = graph.getVertexList();
            for (int v : vert) {
                std::vector<int> edges = graph.getEdges(v);
                std::set<int> l1Wrong;
                std::map<int, std::vector<int> > l2CntWrong;

                for (int e : edges) {
                    int u = graph.getEdgeTo(e);
                    if (aligInfo.isCorrectEdge(&graph, e) == alig_info::InfoAboutContigsAlig::OK) continue;
                    int lib = graph.getEdgeLib(e);

                    if (lib == wlib) {
                        l1Wrong.insert(u);
                    }

                    for (int l : rlibs) {
                        if (l == lib) {
                            l2CntWrong[u].push_back(l);
                        }
                    }
                }

                bool needOut = false;

                for (int u : l1Wrong) {
                    std::sort(l2CntWrong[u].begin(), l2CntWrong[u].end());
                    l2CntWrong[u].resize(std::unique(l2CntWrong[u].begin(), l2CntWrong[u].end())  - l2CntWrong[u].begin());

                    if (l2CntWrong[u].size() < rlibs.size()) {
                        needOut = true;
                    }
                }
                if (needOut) {
                    out << v << " ";
                }
            }

            out.close();
        }
    }
}
