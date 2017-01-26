#ifndef SCAFFOLDER_COMMANDRESETIGNORE_H
#define SCAFFOLDER_COMMANDRESETIGNORE_H

#include <Filter/Filters/Filter.h>
#include "Command.h"

class CommandResetIgnore : public Command {
public:
    void execute(std::string argv, State state, Filter *filter) override;
};


#endif //SCAFFOLDER_COMMANDRESETIGNORE_H
