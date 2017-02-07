#include <iostream>
#include "WriteLocal.h"

WriteLocal::WriteLocal(int v, int dist, std::string fileName, Filter *filter1) : Writer(filter1), v(v),
                                                                                 dist(dist), fileName(fileName) {}

void WriteLocal::write() {
    dotWriter.writeVertexSet(searcher.findVertInLocalArea(v, dist), fileName);
}
