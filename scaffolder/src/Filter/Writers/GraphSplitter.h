#ifndef SCAFFOLDER_GRAPHSPLITTER_H
#define SCAFFOLDER_GRAPHSPLITTER_H

#include <Filter/Filters/Filter.h>
#include <vector>

//split graph on small parts
class GraphSplitter {
private:
    int maxEdge = 15; //max edges in one component
    int maxVert = 10; //max vertexs in one component

    //0 - vertex was not used, 1 - current component, 2 - vertex was used before, 3 - was not used,
    //but it is bed vertex for this component
    std::vector<int> used;
    std::vector< std::vector<int> > res; //res[i] - vertexs in i's component.
    std::vector<int> edgeCol; //component for i-s edge.

    Filter *filter;
    std::vector<int> vert;

    void clear();

    void findNewComp(int v);
public:
    //in graph filter split this set of vertex
    std::vector< std::vector<int> > split(Filter* filter, std::vector<int> vert);
};


#endif //SCAFFOLDER_GRAPHSPLITTER_H
