#ifndef SCAFFOLDER_WRITELOCAL_H
#define SCAFFOLDER_WRITELOCAL_H

#include "Writer.h"

//write only local area of vertex v
class WriteLocal : public  Writer {
private:
    int v; //vertex id, center of local area
    int dist; //dist of vertex to vertex v
    std::string fileName; //file name for writing graph
public:
    WriteLocal(int v, int dist, std::string fileName, Filter *filter1, FileValidator * validator,  int maxVert, int maxEdge);

    void write() override;
};


#endif //SCAFFOLDER_WRITELOCAL_H
