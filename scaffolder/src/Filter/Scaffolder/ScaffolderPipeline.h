#ifndef SCAFFOLDER_SCAFFOLDERPIPELINE_H
#define SCAFFOLDER_SCAFFOLDERPIPELINE_H

#include <Filter/Filters/Filter.h>
#include <string>
#include "Scaffolds.h"

class ScaffolderPipeline {
private:
    Scaffolds scaffolds;

    std::vector<int> topsort;
    std::vector<int> color;
    std::vector<int> topSortPos;

    void topSortDfs(int v, Filter *graph, std::vector<int>* used);
    void colorDfs(int v, int col, Filter * graph);

    bool isUniquePair(int v1, int v2, Filter * graph);
    void topSort(Filter *graph);
    void findCycle(Filter *graph);
    void uniqueConnection(Filter *graph);

    void inOneLineConnection(Filter *graph);
    void addFirstConnection(Filter *graph);
    void delEdgeFromDifPath(Filter *graph);

public:
    void evaluate(Filter* graph, std::string out);

    ScaffolderPipeline(std::string contigFile);

};


#endif //SCAFFOLDER_SCAFFOLDERPIPELINE_H
