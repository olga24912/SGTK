#ifndef SCAFFOLDER_COMMANDHISTOGRAM_H
#define SCAFFOLDER_COMMANDHISTOGRAM_H


#include <Filter/CommandParsers/Commands/Command.h>

class CommandHistogram : public Command {
public:
    void execute(std::string argv, State &state, Filter *filter) override;
};


#endif //SCAFFOLDER_COMMANDHISTOGRAM_H
