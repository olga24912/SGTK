#ifndef SCAFFOLDER_SCAFFOLDSTRATEGYLEAVECONNECTION_H
#define SCAFFOLDER_SCAFFOLDSTRATEGYLEAVECONNECTION_H


#include "ScaffoldStrategy.h"

class ScaffoldStrategyLeaveConnection : public ScaffoldStrategy {
protected:
    void addConnectionForLeaves(Scaffolds *scaffolds, Filter *graph, int dir, std::vector<int> minW);

public:
    void addConnection(Scaffolds *scaffolds, Filter *graph, std::vector<int> minW) override;
};


#endif //SCAFFOLDER_SCAFFOLDSTRATEGYLEAVECONNECTION_H
