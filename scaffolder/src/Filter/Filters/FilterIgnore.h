#ifndef SCAFFOLDER_FILTERIGNORE_H
#define SCAFFOLDER_FILTERIGNORE_H

#include "Filter.h"

class FilterIgnore : public Filter {
private:
    std::vector<bool> ignore;
public:
    FilterIgnore(Filter* filter);

    std::vector<int> getEdges(int v) override;

    std::vector<int> getEdgesR(int v) override;

    std::vector<int> getVertexList() override;

    void processQuery(Query query) override;
};


#endif //SCAFFOLDER_FILTERIGNORE_H
