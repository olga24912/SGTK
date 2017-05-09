#ifndef SCAFFOLDER_COMMANDDIFSTATISTIC_H
#define SCAFFOLDER_COMMANDDIFSTATISTIC_H

#include <Filter/CommandParsers/Commands/Command.h>

class CommandDifStatistic : public Command {
public:
    void execute(std::string argv, State &state, Filter *filter) override;
};


#endif //SCAFFOLDER_COMMANDDIFSTATISTIC_H
