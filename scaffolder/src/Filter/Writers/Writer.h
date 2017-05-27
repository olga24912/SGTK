#ifndef SCAFFOLDER_WRITER_H
#define SCAFFOLDER_WRITER_H

#include <Filter/Filters/Filter.h>
#include <Filter/Writers/DotWriter/DotWriterBuilder.h>
#include "Searcher.h"
#include "Filter/Writers/DotWriter/DotWriter.h"

//abstract class for classes that write graph
class Writer {
protected:
    Filter* filter; //graph for writing
    Searcher searcher = Searcher(filter); //for searching in graph
    DotWriter dotWriter = DotWriter(); //for write graph in dot format
public:
    Writer(Filter* filter1, FileValidator* validator, DotWriterBuilder builder) {
        filter = filter1;
        searcher = Searcher(filter1);
        dotWriter = builder.build();
    }

    virtual void write() = 0;
};


#endif //SCAFFOLDER_WRITER_H
