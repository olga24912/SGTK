//
// Created by olga on 26.01.17.
//

#include "writeLocal.h"

writeLocal::writeLocal(int v, int dist, std::string fileName, Filter *filter1) : Writer(filter1), v(v),
                                                                                 dist(dist), fileName(fileName) {}

void writeLocal::write() {
    dotWriter.writeVertexSet(searcher.findVertInLocalArea(v, dist), fileName);
}
