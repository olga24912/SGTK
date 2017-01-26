#ifndef SCAFFOLDER_SEARCHER_H
#define SCAFFOLDER_SEARCHER_H

#include "Filter/Filters/Filter.h"
#include <algorithm>

class Searcher {
private:
    Filter filter;

    void dfsFindComponent(int v, int currentCol, int* color);
public:
    Searcher(Filter filter1): filter(filter1) {}
    std::vector<int> findVertInLocalArea(int v, int dist);
    int findComponent(int *col);
};


#endif //SCAFFOLDER_SEARCHER_H
