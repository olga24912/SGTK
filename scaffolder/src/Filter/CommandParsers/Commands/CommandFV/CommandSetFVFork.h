#ifndef SCAFFOLDER_COMMANDSETFVFORK_H
#define SCAFFOLDER_COMMANDSETFVFORK_H

#include <Filter/Writers/FileValidator/ValidatorNotPathWithAllLib.h>
#include <Filter/Writers/FileValidator/ValidatorFork.h>
#include "CommandSetFV.h"

class CommandSetFVFork : public CommandSetFV {
protected:
    void setFV(State &state, std::string argv) override {
        std::stringstream ss(argv);
        int ln;
        ss >> ln;
        state.validator = new ValidatorFork(ln);
    }
};


#endif //SCAFFOLDER_COMMANDSETFVFORK_H
