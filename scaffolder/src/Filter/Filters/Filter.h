#ifndef SCAFFOLDER_FILTER_H
#define SCAFFOLDER_FILTER_H

#include <vector>
#include <string>

class Filter {
    Filter* subfilter = nullptr;
    virtual int getVertexCount() = 0;
    virtual std::vector<int> getEdges(int v) = 0;
    virtual std::vector<int> getEdgesR(int v) = 0;
    virtual std::string getTargetName(int v) = 0;
    virtual int getTargetLen(int v) = 0;
    virtual int getEdgeTo(int e) = 0;
    virtual int getEdgeFrom(int e) = 0;
    virtual int getEdgeWieght(int e) = 0;
    virtual std::string getEdgeColor(int e) = 0;
    virtual std::string getEdgeLibName(int e) = 0;
};


#endif //SCAFFOLDER_FILTER_H
