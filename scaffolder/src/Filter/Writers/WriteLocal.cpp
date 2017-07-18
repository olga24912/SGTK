#include <iostream>
#include "WriteLocal.h"

namespace filter {
    namespace writers {
        WriteLocal::WriteLocal(int v, int dist, std::string fileName, ContigGraph *graph1,
                               FileValidator *validator, DotWriterBuilder *builder) :
                Writer(graph1, validator, builder), v(v), dist(dist), fileName(fileName) {}

        void WriteLocal::write() {
            INFO("Write Local");
            dotWriter->writeVertexSet(searcher.findVertInLocalArea(v, dist), fileName);
        }
    }
}