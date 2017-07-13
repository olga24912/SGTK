#ifndef SCAFFOLDER_COMMANDSETFVERROR_H
#define SCAFFOLDER_COMMANDSETFVERROR_H

#include <Filter/Writers/FileValidator/ValidatorError.h>
#include "CommandSetFV.h"

namespace filter {
    namespace commands {
        class CommandSetFVError : public CommandSetFV {
        protected:
            void setFV(State &state, std::string argv) override {
                INFO("set file validator Error");
                std::stringstream ss(argv);
                std::string fileName;
                int libNum;
                ss >> fileName >> libNum;

                state.validator = new writers::ValidatorError(fileName, libNum);
            }
        };
    }
}

#endif //SCAFFOLDER_COMMANDSETFVERROR_H
