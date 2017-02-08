#ifndef SCAFFOLDER_SIMPLEWRITER_H
#define SCAFFOLDER_SIMPLEWRITER_H


#include <Filter/Filters/Filter.h>
#include <fstream>
#include <algorithm>
#include "GraphSplitter.h"

//class for write graph in dot format.
class DotWriter {
private:
    Filter* filter; //graph for writing
    GraphSplitter graphSplitter; //split graph on small parts

    void writeOneVertex(int v, bool isColored, std::ofstream& out);
    void writeOneEdge(int e, std::ofstream& out);

    void writeOneVertexSet(std::vector<int> vert, std::string fileName);

    bool isGoodFileForWrite(std::vector<int> vert);
public:
    DotWriter(Filter* filter): filter(filter){}

    //write this set of vertex in files with prefix fileName
    void writeVertexSet(std::vector<int> vert, std::string fileName);
};


#endif //SCAFFOLDER_SIMPLEWRITER_H
