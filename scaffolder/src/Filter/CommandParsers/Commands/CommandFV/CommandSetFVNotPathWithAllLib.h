#ifndef SCAFFOLDER_COMMANDSETFVNOTPATHWITHALLLIB_H
#define SCAFFOLDER_COMMANDSETFVNOTPATHWITHALLLIB_H

#include <Filter/Writers/FileValidator/ValidatorNotPathWithAllLib.h>
#include "CommandSetFV.h"

class CommandSetFVNotWithAllLib : public CommandSetFV {
protected:
    void setFV(State &state, std::string argv) override {
        state.validator = new ValidatorNotPathWithAllLib();
    }
};

#endif //SCAFFOLDER_COMMANDSETFVNOTPATHWITHALLLIB_H
