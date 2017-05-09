#ifndef SCAFFOLDER_COMMANDWEIGHTSTATISTIC_H
#define SCAFFOLDER_COMMANDWEIGHTSTATISTIC_H

#include <Filter/Statistics/Statistic.h>
#include <Filter/CommandParsers/State.h>
#include <Filter/CommandParsers/Commands/Command.h>

class CommandWeightStatistic : public Command {
public:
    void execute(std::string argv, State &state, Filter *filter) override;
};


#endif //SCAFFOLDER_COMMANDWEIGHTSTATISTIC_H
