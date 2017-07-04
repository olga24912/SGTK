#ifndef SCAFFOLDER_COMMANDSETFVONLYFIRST_H
#define SCAFFOLDER_COMMANDSETFVONLYFIRST_H

#include <Filter/Writers/FileValidator/ValidatorNotPathWithAllLib.h>
#include <Filter/Writers/FileValidator/ValidatorFork.h>
#include <Filter/Writers/FileValidator/ValidatorOnlyFirst.h>
#include "CommandSetFV.h"

class CommandSetFVOnlyFirst : public CommandSetFV {
protected:
    void setFV(State &state, std::string argv) override {
        std::stringstream ss(argv);
        int ln1, ln2;
        ss >> ln1 >> ln2;

        INFO("set file validator with first" << ln1 << " but without second " << ln2 << " lib");

        state.validator = new ValidatorOnlyFirst(ln1, ln2);
    }
};

#endif //SCAFFOLDER_COMMANDSETFVONLYFIRST_H
