#ifndef SCAFFOLDER_FILTERIGNOREEDGE_H
#define SCAFFOLDER_FILTERIGNOREEDGE_H


#include <vector>
#include "Filter.h"

class FilterIgnoreEdge : public Filter {
private:
    std::vector<bool> ignore;
public:
    FilterIgnoreEdge(Filter* filter);

    std::vector<int> getEdges(int v) override;

    std::vector<int> getEdgesR(int v) override;

    //handling querys SET_IGNORE_EDGE with args <edge id for ignore>
    void processQuery(Query query) override;
};

#endif //SCAFFOLDER_FILTERIGNOREEDGE_H
