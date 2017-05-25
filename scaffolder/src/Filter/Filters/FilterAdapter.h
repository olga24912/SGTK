#ifndef SCAFFOLDER_FILTERADAPTER_H
#define SCAFFOLDER_FILTERADAPTER_H

#include "Filter.h"
#include "ContigGraph/ContigGraph.h"

//change graph interface to Filter interface
class FilterAdapter: public Filter {
private:
    ContigGraph graph; //graph for change interface
public:
    FilterAdapter(ContigGraph graph): graph(graph) {}

    std::vector<int> getLibList() override;

    int getVertexCount() override;

    std::vector<int> getEdges(int v) override;

    std::vector<int> getEdgesR(int v) override;

    std::string getTargetName(int v) override;

    int getTargetLen(int v) override;

    int getEdgeTo(int e) override;

    int getEdgeFrom(int e) override;

    int getEdgeWeight(int e) override;

    int getEdgeLib(int e) override;

    std::string getLibName(int l) override;

    std::string getLibColor(int l) override;

    std::vector<int> getVertexList() override;
    
    std::string getInfo(int e) override;

    std::pair<int, int> getFirstCoord(int e) override;

    std::pair<int, int> getSecondCoord(int e) override;

    //handling query type UPLOAD_GRAPH with args <contigGraph file name>
    void processQuery(Query query) override;
};


#endif //SCAFFOLDER_FILTERADAPTER_H
