#ifndef SCAFFOLDER_SCAFFOLDSTRATEGYONELINE_H
#define SCAFFOLDER_SCAFFOLDSTRATEGYONELINE_H


#include "ScaffoldStrategy.h"

class ScaffoldStrategyOneLine : public ScaffoldStrategy {
protected:
    void addFirstConnection(Scaffolds *scaffolds, Filter *graph, int minW);
    void delEdgeFromDifPath(Scaffolds *scaffolds, Filter *graph);
public:
    void addConnection(Scaffolds *scaffolds, Filter *graph, int minW) override;

};


#endif //SCAFFOLDER_SCAFFOLDSTRATEGYONELINE_H
