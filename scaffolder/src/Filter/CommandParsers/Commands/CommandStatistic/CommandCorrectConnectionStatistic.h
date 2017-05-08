#ifndef SCAFFOLDER_COMMANDCORRECTCONNECTIONSTATISTIC_H
#define SCAFFOLDER_COMMANDCORRECTCONNECTIONSTATISTIC_H


#include <Filter/CommandParsers/Commands/Command.h>

class CommandCorrectConnectionStatistic : public Command {
public:
    void execute(std::string argv, State &state, Filter *filter) override;
};


#endif //SCAFFOLDER_COMMANDCORRECTCONNECTIONSTATISTIC_H
