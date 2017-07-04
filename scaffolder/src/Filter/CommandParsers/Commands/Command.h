#ifndef SCAFFOLDER_COMMAND_H
#define SCAFFOLDER_COMMAND_H

#include <string>
#include "Filter/CommandParsers/State.h"
#include <sstream>

class Command {
public:
    virtual void execute(std::string argv, State& state, Filter* filter) = 0;
protected:
    DECL_LOGGER("Command");
};


#endif //SCAFFOLDER_COMMAND_H
