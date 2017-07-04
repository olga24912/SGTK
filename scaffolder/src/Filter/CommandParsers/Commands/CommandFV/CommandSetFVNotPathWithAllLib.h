#ifndef SCAFFOLDER_COMMANDSETFVNOTPATHWITHALLLIB_H
#define SCAFFOLDER_COMMANDSETFVNOTPATHWITHALLLIB_H

#include <Filter/Writers/FileValidator/ValidatorNotPathWithAllLib.h>
#include "CommandSetFV.h"

class CommandSetFVNotWithAllLib : public CommandSetFV {
protected:
    void setFV(State &state, std::string argv) override {
        INFO("set file validator not path with all libs");
        state.validator = new ValidatorNotPathWithAllLib();
    }
};

#endif //SCAFFOLDER_COMMANDSETFVNOTPATHWITHALLLIB_H
