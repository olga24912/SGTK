#ifndef SCAFFOLDER_SIMPLEWRITER_H
#define SCAFFOLDER_SIMPLEWRITER_H


#include <Filter/Filters/Filter.h>
#include <fstream>
#include <algorithm>
#include <Filter/Writers/FileValidator/FileValidator.h>
#include <Filter/Writers/FileValidator/ValidatorNotPathWithAllLib.h>
#include "Filter/Writers/GraphSplitter.h"

//class for write graph in dot format.
class DotWriter {
protected:
    Filter* filter; //graph for writing
    GraphSplitter graphSplitter; //split graph on small parts
    FileValidator* validator = new ValidatorNotPathWithAllLib();

    virtual void writeOneVertex(int v, bool isColored, std::ofstream& out);

    virtual void writeOneEdge(int e, std::ofstream& out);

    virtual void writeOneVertexSet(std::vector<int> vert, std::string fileName);

public:
    DotWriter(){}
    DotWriter(Filter* filter): filter(filter){}
    DotWriter(Filter* filter, FileValidator* validator, int maxVert, int maxEdge):
            filter(filter), validator(validator){
        graphSplitter = GraphSplitter(maxVert, maxEdge);
    }

    //write this set of vertex in files with prefix fileName
    virtual void writeVertexSet(std::vector<int> vert, std::string fileName);

protected:
    DECL_LOGGER("DotWriter");
};


#endif //SCAFFOLDER_SIMPLEWRITER_H
