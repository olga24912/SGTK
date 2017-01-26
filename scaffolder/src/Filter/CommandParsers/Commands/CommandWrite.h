#ifndef SCAFFOLDER_COMMANDWRITE_H
#define SCAFFOLDER_COMMANDWRITE_H

#include <Filter/Filters/Filter.h>
#include <Filter/CommandParsers/State.h>
#include <Filter/Tools/SystemTools.h>
#include "Command.h"

class CommandWrite : public Command {
public:
    void execute(std::string argv, State& state, Filter *filter) override;

    virtual void writeGraph(std::string argv, State& state, Filter *filter) = 0;
};


#endif //SCAFFOLDER_COMMANDWRITE_H
