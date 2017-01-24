#ifndef SCAFFOLDER_FILTER_H
#define SCAFFOLDER_FILTER_H

#include <vector>
#include <string>
#include "assert.h"

class Filter {
    Filter* subfilter = nullptr;
    virtual int getVertexCount();
    virtual std::vector<int> getEdges(int v);
    virtual std::vector<int> getEdgesR(int v);
    virtual std::string getTargetName(int v);
    virtual int getTargetLen(int v);
    virtual int getEdgeTo(int e);
    virtual int getEdgeFrom(int e);
    virtual int getEdgeWieght(int e);
    virtual std::string getEdgeColor(int e);
    virtual std::string getEdgeLibName(int e);
};


#endif //SCAFFOLDER_FILTER_H
