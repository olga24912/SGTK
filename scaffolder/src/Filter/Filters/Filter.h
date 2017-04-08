#ifndef SCAFFOLDER_FILTER_H
#define SCAFFOLDER_FILTER_H

#include <vector>
#include <string>
#include "assert.h"
#include "Query.h"
#include <sstream>

//abstract class for all filters
class Filter {
protected:
    Filter(){};
    Filter(Filter *filer): subfilter(filer) {}
    Filter* subfilter = nullptr; //graph for filtering
public:
    virtual std::vector<int> getLibList(); //return list of libs num with filters
    virtual std::vector<int> getVertexList(); //return filter list of vertexs id
    virtual int getVertexCount(); //return count of vertexs without filters
    virtual std::vector<int> getEdges(int v); //return filter edges id from v
    virtual std::vector<int> getEdgesR(int v); //return filter edges id to v
    virtual std::string getTargetName(int v); //return contig name of vertex v
    virtual int getTargetLen(int v); //return contig len of vertex v
    virtual int getEdgeTo(int e); //if e is edge v->u, return u
    virtual int getEdgeFrom(int e); //if e is edge v->u, return v
    virtual int getEdgeWeight(int e); //return edge weight for edge with id e
    virtual int getEdgeLib(int e); //return edge lib num
    virtual std::string getLibName(int l); //return name for lib
    virtual std::string getLibColor(int l); //return color for this lib

    virtual void processQuery(Query query); //handling query for changing filter options

    void write(std::string fileName); //serialize this graph in .gr format in "fileName" file

    ~Filter() {
        if (subfilter != nullptr) {
            delete subfilter;
        }
    }
};


#endif //SCAFFOLDER_FILTER_H
