#ifndef SCAFFOLDER_WRITER_H
#define SCAFFOLDER_WRITER_H

#include <Filter/Filters/ContigGraph.h>
#include <Filter/Writers/DotWriter/DotWriterBuilder.h>
#include "Searcher.h"
#include "Filter/Writers/DotWriter/DotWriter.h"

namespace filter {
    namespace writers {
        using namespace contig_graph;
//abstract class for classes that write graph
        class Writer {
        protected:
            ContigGraph *graph; //graph for writing
            Searcher searcher = Searcher(graph); //for searching in graph
            DotWriter *dotWriter; //for write graph in dot format
        public:
            Writer(ContigGraph *graph1, FileValidator *validator, DotWriterBuilder *builder) :
                    graph(graph1), searcher(Searcher(graph1)), dotWriter(builder->build()) {
            }

            virtual void write() = 0;

            ~Writer() {
                delete dotWriter;
            }
        };
    }
}

#endif //SCAFFOLDER_WRITER_H
