#ifndef SCAFFOLDER_COMMANDFVWITHDIFINLIB_H
#define SCAFFOLDER_COMMANDFVWITHDIFINLIB_H

#include "CommandSetFV.h"

class CommandFVWithDifInLib : public CommandSetFV {
protected:
    void setFV(State &state, std::string argv) override;
};


#endif //SCAFFOLDER_COMMANDFVWITHDIFINLIB_H
