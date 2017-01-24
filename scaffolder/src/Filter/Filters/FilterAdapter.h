#ifndef SCAFFOLDER_FILTERADAPTER_H
#define SCAFFOLDER_FILTERADAPTER_H

#include "Filter.h"
#include "ContigGraph/ContigGraph.h"

class FilterAdapter: public Filter {
    ContigGraph graph;

    FilterAdapter(ContigGraph graph): graph(graph) {}

    int getVertexCount() override;

    std::vector<int> getEdges(int v) override;

    std::vector<int> getEdgesR(int v) override;

    std::string getTargetName(int v) override;

    int getTargetLen(int v) override;

    int getEdgeTo(int e) override;

    int getEdgeFrom(int e) override;

    int getEdgeWieght(int e) override;

    std::string getEdgeColor(int e) override;

    std::string getEdgeLibName(int e) override;
};


#endif //SCAFFOLDER_FILTERADAPTER_H
