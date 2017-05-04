#ifndef SCAFFOLDER_SCAFFOLDERPIPELINE_H
#define SCAFFOLDER_SCAFFOLDERPIPELINE_H

#include <Filter/Filters/Filter.h>
#include <string>
#include "Scaffolds.h"

class ScaffolderPipeline {
public:
    void evaluate(Filter* graph, std::string contigFile, std::string out);
};


#endif //SCAFFOLDER_SCAFFOLDERPIPELINE_H
