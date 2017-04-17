#ifndef SCAFFOLDER_COMMANDMINCONTIG_H
#define SCAFFOLDER_COMMANDMINCONTIG_H

#include <Filter/Filters/Filter.h>
#include "Filter/CommandParsers/Commands/Command.h"

class CommandMinContig : public Command {
public:
    void execute(std::string argv, State& state, Filter *filter) override;
};


#endif //SCAFFOLDER_COMMANDMINCONTIG_H
