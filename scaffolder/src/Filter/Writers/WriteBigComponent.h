#ifndef SCAFFOLDER_WRITEBIGCOMPONENT_H
#define SCAFFOLDER_WRITEBIGCOMPONENT_H

#include "Writer.h"

class WriteBigComponent : public Writer {
private:
    std::string fileName;
    int minSize;
public:
    WriteBigComponent(std::string fileName, int minSize, Filter *filter1);

    void write() override;

};


#endif //SCAFFOLDER_WRITEBIGCOMPONENT_H
