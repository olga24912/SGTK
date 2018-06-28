#ifndef SCAFFOLDER_COMMANDSETFVNOTSIMPLEPATH_H
#define SCAFFOLDER_COMMANDSETFVNOTSIMPLEPATH_H

#include <Filter/Writers/FileValidator/ValidatorNotSimplePath.h>
#include "CommandSetFV.h"

namespace filter {
    namespace commands {
        class CommandSetFVNotSimplePath : public CommandSetFV {
        protected:
            void setFV(State &state, std::string argv) override {
                INFO("set file validator");
                state.validator = new writers::ValidatorNotSimplePath();
            }
        };
    }
}

#endif //SCAFFOLDER_COMMANDSETFVNOTSIMPLEPATH_H
