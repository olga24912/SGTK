#ifndef SCAFFOLDER_COMMANDMERGESIMPLEPATH_H
#define SCAFFOLDER_COMMANDMERGESIMPLEPATH_H

#include <Filter/Filters/Filter.h>
#include <Filter/CommandParsers/State.h>
#include "Filter/CommandParsers/Commands/Command.h"

class CommandScaffold : public Command {
public:
    void execute(std::string argv, State &state, Filter *filter) override;
};


#endif //SCAFFOLDER_COMMANDMERGESIMPLEPATH_H
