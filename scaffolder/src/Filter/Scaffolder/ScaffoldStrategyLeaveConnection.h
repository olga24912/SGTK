#ifndef SCAFFOLDER_SCAFFOLDSTRATEGYLEAVECONNECTION_H
#define SCAFFOLDER_SCAFFOLDSTRATEGYLEAVECONNECTION_H


#include "ScaffoldStrategy.h"

class ScaffoldStrategyLeaveConnection : public ScaffoldStrategy {
protected:
    void addConnectionForLeaves(Scaffolds *scaffolds, Filter *graph, int dir, int minW);

public:
    void addConnection(Scaffolds *scaffolds, Filter *graph, int minW) override;
};


#endif //SCAFFOLDER_SCAFFOLDSTRATEGYLEAVECONNECTION_H
