#include <algorithm>
#include <iostream>
#include <set>
#include "ScaffolderPipeline.h"
#include "ScaffoldStrategy.h"
#include "ScaffoldStrategyOneLine.h"
#include "ScaffoldStrategyLeaveConnection.h"

void ScaffolderPipeline::evaluate(Filter *graph, std::string contigFile, std::string out) {
    Scaffolds scaffolds(contigFile);

    std::vector<int> libs = graph->getLibList();

    std::vector<int> minEdge((unsigned long) libs[(int)libs.size() - 1], 1e9);
    std::vector<int> vert = graph->getVertexList();
    for (int v : vert) {
        std::vector<int> edges = graph->getEdges(v);
        for (int e : edges) {
            minEdge[graph->getEdgeLib(e)] = std::min(minEdge[graph->getEdgeLib(e)], graph->getEdgeWeight(e));
        }
    }

   for (int w = 1; w < 8; w *= 2) {
        ScaffoldStrategyOneLine scafol;
        ScaffoldStrategyLeaveConnection scafl;


        for (int l : libs) {
            std::stringstream arg;
            arg << l << " " << std::max(w, minEdge[l]);
            std::cerr << arg.str() << std::endl;
            graph->processQuery(Query(Query::MIN_EDGE_WEIGHT, arg.str()));
        }

        scafol.addConnection(&scaffolds, graph);
        scafl.addConnection(&scaffolds, graph);

        std::cerr << w << " " << scaffolds.lineId(2501) << " " << scaffolds.isFirst(2501) << " " << scaffolds.isLast(2501) << std::endl;
   }

    scaffolds.print(out);
}
