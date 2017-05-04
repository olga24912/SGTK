#ifndef SCAFFOLDER_SCAFFOLDSTRATEGYUNIQUECONNECTION_H
#define SCAFFOLDER_SCAFFOLDSTRATEGYUNIQUECONNECTION_H


#include "ScaffoldStrategy.h"

class ScaffoldStrategyUniqueConnection : public ScaffoldStrategy {
private:
    bool isUniquePair(int v1, int v2, Filter * graph);
public:
    void addConnection(Scaffolds *scaffolds, Filter *graph) override;
};


#endif //SCAFFOLDER_SCAFFOLDSTRATEGYUNIQUECONNECTION_H
