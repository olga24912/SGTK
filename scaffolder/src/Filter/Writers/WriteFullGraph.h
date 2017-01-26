#ifndef SCAFFOLDER_WRITEFULLGRAPH_H
#define SCAFFOLDER_WRITEFULLGRAPH_H

#include "Writer.h"

class WriteFullGraph : public Writer {
private:
    std::string fileName;
public:
    WriteFullGraph(std::string fileName, Filter *filter1);

    void write() override;
};


#endif //SCAFFOLDER_WRITEFULLGRAPH_H
