#ifndef SCAFFOLDER_WRITER_H
#define SCAFFOLDER_WRITER_H

#include <Filter/Filters/Filter.h>
#include "Searcher.h"
#include "DotWriter.h"

//abstract class for classes that write graph
class Writer {
protected:
    Filter* filter; //graph for writing
    Searcher searcher = Searcher(filter); //for searching in graph
    DotWriter dotWriter = DotWriter(filter); //for write graph in dot format
public:
    Writer(Filter* filter1, FileValidator* validator) {
        filter = filter1;
        searcher = Searcher(filter1);
        dotWriter = DotWriter(filter1, validator);
    }

    virtual void write() = 0;
};


#endif //SCAFFOLDER_WRITER_H
