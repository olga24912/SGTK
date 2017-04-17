#ifndef SCAFFOLDER_WRITEFULLGRAPH_H
#define SCAFFOLDER_WRITEFULLGRAPH_H

#include "Writer.h"

//Write full graph
class WriteFullGraph : public Writer {
private:
    std::string fileName; //fileName for write full graph
public:
    WriteFullGraph(std::string fileName, Filter *filter1, FileValidator *validator,
                   int maxVert, int maxEdge);

    void write() override;
};


#endif //SCAFFOLDER_WRITEFULLGRAPH_H
