#ifndef SCAFFOLDER_SIMPLEWRITER_H
#define SCAFFOLDER_SIMPLEWRITER_H


#include <Filter/Filters/Filter.h>
#include <fstream>
#include <algorithm>
#include "GraphSplitter.h"

class DotWriter {
private:
    Filter* filter;
    GraphSplitter graphSplitter;

    void writeOneVertex(int v, bool isColored, std::ofstream& out);
    void writeOneEdge(int e, std::ofstream& out);

    void writeOneVertexSet(std::vector<int> vert, std::string fileName);
public:
    DotWriter(Filter* filter): filter(filter){}

    void writeVertexSet(std::vector<int> vert, std::string fileName);

};


#endif //SCAFFOLDER_SIMPLEWRITER_H
