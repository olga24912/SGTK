#ifndef SCAFFOLDER_COMMANDMINEDGEWEIGHT_H
#define SCAFFOLDER_COMMANDMINEDGEWEIGHT_H

#include <Filter/Filters/Filter.h>
#include "Command.h"

class CommandMinEdgeWeight : public Command {
public:
    void execute(std::string argv, State& state, Filter *filter) override;
};


#endif //SCAFFOLDER_COMMANDMINEDGEWEIGHT_H
