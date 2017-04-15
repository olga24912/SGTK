#ifndef SCAFFOLDER_FILEVALIDATOR_H
#define SCAFFOLDER_FILEVALIDATOR_H


#include <vector>
#include <Filter/Filters/Filter.h>

class FileValidator {
public:
    virtual bool isGoodVertexSet(std::vector<int> vert, Filter* filter) {
        return true;
    }
};


#endif //SCAFFOLDER_FILEVALIDATOR_H
