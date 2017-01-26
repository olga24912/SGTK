#ifndef SCAFFOLDER_WRITER_H
#define SCAFFOLDER_WRITER_H

#include <Filter/Filters/Filter.h>
#include "Searcher.h"
#include "DotWriter.h"

class Writer {
protected:
    Filter* filter;
    Searcher searcher = Searcher(filter);
    DotWriter dotWriter = DotWriter(filter);
public:
    Writer(Filter* filter1) {
        filter = filter1;
        searcher = Searcher(filter1);
        dotWriter = DotWriter(filter1);
    }

    virtual void write() = 0;
};


#endif //SCAFFOLDER_WRITER_H
