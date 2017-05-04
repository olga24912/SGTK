#ifndef SCAFFOLDER_SCAFFOLDSTRATEGY_H
#define SCAFFOLDER_SCAFFOLDSTRATEGY_H


#include <vector>
#include <Filter/Filters/Filter.h>
#include "Scaffolds.h"

class ScaffoldStrategy {
protected:
    std::vector<int> topsort;
    std::vector<int> color;
    std::vector<int> topSortPos;

    void topSortDfs(int v, Filter *graph, std::vector<int>* used);
    void colorDfs(int v, int col, Filter * graph);

    void topSort(Filter *graph);
    void findCycle(Filter *graph);

    int deg(int i, Filter *pFilter, int dirIn);

public:
    virtual void addConnection(Scaffolds* scaffolds, Filter *graph) = 0;
};


#endif //SCAFFOLDER_SCAFFOLDSTRATEGY_H
