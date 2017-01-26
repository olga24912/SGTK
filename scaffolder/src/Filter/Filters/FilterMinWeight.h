#ifndef SCAFFOLDER_FILTERMINWEIDHT_H
#define SCAFFOLDER_FILTERMINWEIDHT_H

#include "Filter.h"

class FilterMinWeight: public Filter {
private:
    int minContigLen = 500;
    std::vector<int> libMinEdgeWight;
public:
    FilterMinWeight(Filter* filter);

    std::vector<int> getEdges(int v) override;

    std::vector<int> getEdgesR(int v) override;

    std::vector<int> getVertexList() override;

    void processQuery(Query query) override;
};


#endif //SCAFFOLDER_FILTERMINWEIDHT_H
