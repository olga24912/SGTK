#ifndef SCAFFOLDER_FILTER_H
#define SCAFFOLDER_FILTER_H

#include <vector>
#include <string>
#include "assert.h"
#include "Query.h"
#include <sstream>

class Filter {
protected:
    Filter(){};
    Filter(Filter *filer): subfilter(filer) {}
    Filter* subfilter = nullptr;
public:
    virtual int getLibCount();
    virtual std::vector<int> getVertexList();
    virtual int getVertexCount();
    virtual std::vector<int> getEdges(int v);
    virtual std::vector<int> getEdgesR(int v);
    virtual std::string getTargetName(int v);
    virtual int getTargetLen(int v);
    virtual int getEdgeTo(int e);
    virtual int getEdgeFrom(int e);
    virtual int getEdgeWieght(int e);
    virtual int getEdgeLib(int e);
    virtual std::string getLibName(int l);
    virtual std::string getLibColor(int l);

    virtual void processQuery(Query query);

    ~Filter() {
        if (subfilter != nullptr) {
            delete subfilter;
        }
    }
};


#endif //SCAFFOLDER_FILTER_H
