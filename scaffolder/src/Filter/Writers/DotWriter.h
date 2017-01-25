//
// Created by olga on 25.01.17.
//

#ifndef SCAFFOLDER_SIMPLEWRITER_H
#define SCAFFOLDER_SIMPLEWRITER_H


#include <Filter/Filters/Filter.h>
#include <fstream>
#include <algorithm>

class DotWriter {
private:
    Filter filter;

    void writeOneVertex(int v, bool isColored, std::ofstream& out);
    void writeOneEdge(int e, std::ofstream& out);
public:
    DotWriter(Filter filter): filter(filter){}

    void writeVertexSet(std::vector<int> vert, std::string fileName);

};


#endif //SCAFFOLDER_SIMPLEWRITER_H
