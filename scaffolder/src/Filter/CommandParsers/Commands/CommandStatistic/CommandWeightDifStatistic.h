#ifndef SCAFFOLDER_COMMANDWEIGHTDIFSTATISTIC_H
#define SCAFFOLDER_COMMANDWEIGHTDIFSTATISTIC_H

#include <Filter/CommandParsers/Commands/Command.h>

class CommandWeightDifStatistic : public Command {
public:
    void execute(std::string argv, State &state, Filter *filter) override;
};


#endif //SCAFFOLDER_COMMANDWEIGHTDIFSTATISTIC_H
