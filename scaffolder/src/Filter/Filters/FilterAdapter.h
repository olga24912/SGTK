#ifndef SCAFFOLDER_FILTERADAPTER_H
#define SCAFFOLDER_FILTERADAPTER_H

#include "Filter.h"
#include "ContigGraph/ContigGraph.h"

class FilterAdapter: public Filter {
private:
    ContigGraph graph;
public:
    FilterAdapter(ContigGraph graph): graph(graph) {}

    int getLibCount() override;

    int getVertexCount() override;

    std::vector<int> getEdges(int v) override;

    std::vector<int> getEdgesR(int v) override;

    std::string getTargetName(int v) override;

    int getTargetLen(int v) override;

    int getEdgeTo(int e) override;

    int getEdgeFrom(int e) override;

    int getEdgeWieght(int e) override;

    int getEdgeLib(int e) override;

    std::string getLibName(int l) override;

    std::string getLibColor(int l) override;

    std::vector<int> getVertexList() override;

    void processQuery(Query query) override;
};


#endif //SCAFFOLDER_FILTERADAPTER_H
