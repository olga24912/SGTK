#ifndef SCAFFOLDER_WRITEALONGPATH_H
#define SCAFFOLDER_WRITEALONGPATH_H

#include "Writer.h"

//write only edges and vertex in local area from some path
class WriteAlongPath : public Writer {
private:
    std::string fileName; //file for write the graph in dot format
    int libId; //lib id for find path using only this edges
    int dist; //area from path
    int minSize; //minLen of path for writing
public:
    WriteAlongPath(std::string fileName, int libId, int dist, int minSize, Filter *filter1, FileValidator *validator,
                   int maxVert, int maxEdge);

    void write() override;
};


#endif //SCAFFOLDER_WRITEALONGPATH_H
