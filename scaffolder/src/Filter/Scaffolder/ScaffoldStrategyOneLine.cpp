#include "ScaffoldStrategyOneLine.h"

void ScaffoldStrategyOneLine::addConnection(Scaffolds *scaffolds, Filter *graph, std::vector<int> minW) {
    newCon.resize(graph->getVertexCount(), 0);
    newConR.resize(graph->getVertexCount(), 0);
    topSort(graph);
    findCycle(graph);
    addFirstConnection(scaffolds, graph, minW);
    delEdgeFromDifPath(scaffolds, graph);
}

void ScaffoldStrategyOneLine::addFirstConnection(Scaffolds *scaffolds, Filter *graph, std::vector<int> minW) {
    for (int i = 0; i < (int)topsort.size(); ++i) {
        int v = topsort[i];

        std::vector<int> edges = graph->getEdges(v);
        if (edges.size() == 0) continue;
        int minu = -1;

        for (int e : edges) {
            int u = graph->getEdgeTo(e);

            if ((minu == -1 || topSortPos[u] < topSortPos[minu]) && color[v] != color[u] &&
                    graph->getEdgeWeight(e) >= minW[graph->getLibType(graph->getEdgeLib(e))]) {
                minu = u;
            }
        }
        if (v == 10952) {
            std::cerr <<"Min u for 10952: " << color[v] << " " << color[11523] << " " << minu << " " << minW[ContigGraph::Lib::Type::RNA_PAIR] << " " << minW[ContigGraph::Lib::Type::RNA_SPLIT_50] << std::endl;
        }

        if (minu != -1 && color[v] != color[minu] && scaffolds->lineId(v) != scaffolds->lineId(minu) &&
                scaffolds->isLast(v) && scaffolds->isFirst(minu)) {
            scaffolds->addConnection(v, minu);
            if (v == 10952) {
                std::cerr << "add con 10952" << std::endl;
            }
            newCon[v] = 1;
            newConR[minu] = 1;
        }
    }
}

void ScaffoldStrategyOneLine::delEdgeFromDifPath(Scaffolds *scaffolds, Filter *graph) {
    for (int i = (int)topsort.size() - 1; i >= 0; --i) {
        int v = topsort[i];
        std::vector<int> edges = graph->getEdges(v);
        if (edges.size() == 0) continue;

        for (int e : edges) {
            int u = graph->getEdgeTo(e);
            assert(graph->getEdgeFrom(e) == v);

            if(v == 15533) {
                std::cerr << v << " " << u << " " << scaffolds->lineId(v) << " " << scaffolds->lineId(u) <<" " << newCon[v] << std::endl;
            }

            if (scaffolds->lineId(v) != scaffolds->lineId(u) && newCon[v] == 1) {
                if (v == 15533) {
                    std::cerr << "brok con" << v << " " << u << std::endl;
                }

                scaffolds->brokeConnection(v);
            }
            if (scaffolds->lineId(v) != scaffolds->lineId(u) && newConR[u] == 1) {
                if (v == 15533) {
                    std::cerr << "brok con" << v << " " << u << std::endl;
                }

                scaffolds->brokeConnectionTo(u);
            }

        }
    }
}