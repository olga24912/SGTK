#ifndef SCAFFOLDER_COMMANDSETFVFORK_H
#define SCAFFOLDER_COMMANDSETFVFORK_H

#include <Filter/Writers/FileValidator/ValidatorNotPathWithAllLib.h>
#include <Filter/Writers/FileValidator/ValidatorFork.h>
#include "CommandSetFV.h"

namespace filter {
    namespace commands {

        class CommandSetFVFork : public CommandSetFV {
        protected:
            void setFV(State &state, std::string argv) override {
                INFO("set file validator in lib=" << argv);
                std::stringstream ss(argv);
                int ln;
                ss >> ln;
                state.validator = new writers::ValidatorFork(ln);
            }
        };
    }
}

#endif //SCAFFOLDER_COMMANDSETFVFORK_H