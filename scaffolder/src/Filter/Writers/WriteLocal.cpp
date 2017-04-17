#include <iostream>
#include "WriteLocal.h"

WriteLocal::WriteLocal(int v, int dist, std::string fileName, Filter *filter1,
                       FileValidator *validator,  int maxVert, int maxEdge) :
        Writer(filter1, validator, maxVert, maxEdge), v(v), dist(dist), fileName(fileName) {}

void WriteLocal::write() {
    dotWriter.writeVertexSet(searcher.findVertInLocalArea(v, dist), fileName);
}
