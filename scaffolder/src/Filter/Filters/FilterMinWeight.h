#ifndef SCAFFOLDER_FILTERMINWEIDHT_H
#define SCAFFOLDER_FILTERMINWEIDHT_H

#include "Filter.h"

//ignore vertex with small contig len or
//ignore edge with small wieght
class FilterMinWeight: public Filter {
private:
    int minContigLen = 500;
    std::vector<int> libMinEdgeWight;
public:
    FilterMinWeight(Filter* filter);

    std::vector<int> getEdges(int v) override;

    std::vector<int> getEdgesR(int v) override;

    std::vector<int> getVertexList() override;

    //handling query with type MIN_EDGE_WEIGHT with args <lib num> <min edge weight for this lib>
    // and MIN_CONTIG_LEN with args <min contig weight>
    void processQuery(Query query) override;
};


#endif //SCAFFOLDER_FILTERMINWEIDHT_H
