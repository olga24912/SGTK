#ifndef SCAFFOLDER_COMMANDTWOCOMPSTATISTIC_H
#define SCAFFOLDER_COMMANDTWOCOMPSTATISTIC_H

#include <Filter/CommandParsers/Commands/Command.h>

class CommandTwoCompStatistic : public Command {
public:
    void execute(std::string argv, State &state, Filter *filter) override;
};


#endif //SCAFFOLDER_COMMANDTWOCOMPSTATISTIC_H