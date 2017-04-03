#ifndef SCAFFOLDER_COMMANDPRINT_H
#define SCAFFOLDER_COMMANDPRINT_H

#include <Filter/Filters/Filter.h>
#include "Command.h"

class CommandPrint: public Command{

public:
    void execute(std::string argv, State &state, Filter *filter) override;
};


#endif //SCAFFOLDER_COMMANDPRINT_H
