#include <Filter/AligInfo/InfoAboutContigsAlig.h>
#include "RuleDelWrongCoonection.h"

using namespace filter::contig_graph;

namespace filter {
    namespace scaffolder {
        void RuleDelWrongConnection::simplifyGraph(ContigGraph *graph) {
            alig_info::InfoAboutContigsAlig aligInfo;
            aligInfo.parseCoordFile(graph, coordFile);

            std::vector<int> vertList = graph->getVertexList();
            for (int v : vertList) {
                std::vector<int> edges = graph->getEdges(v);
                for (int e : edges) {
                    if (aligInfo.isCorrectEdge(graph, e) != alig_info::InfoAboutContigsAlig::OK) {
                        graph->delEdge(e);
                    }
                }

                edges = graph->getEdgesR(v);
                for (int e : edges) {
                    if (aligInfo.isCorrectEdge(graph, e) != alig_info::InfoAboutContigsAlig::OK) {
                        graph->delEdge(e);
                    }
                }
            }
        }
    }
}
