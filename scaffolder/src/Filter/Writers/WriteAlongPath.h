#ifndef SCAFFOLDER_WRITEALONGPATH_H
#define SCAFFOLDER_WRITEALONGPATH_H

#include "Writer.h"

class WriteAlongPath : public Writer {
private:
    std::string fileName;
    int libId;
    int dist;
    int minSize;
public:
    WriteAlongPath(std::string fileName, int libId, int dist, int minSize, Filter *filter1);

    void write() override;
};


#endif //SCAFFOLDER_WRITEALONGPATH_H
