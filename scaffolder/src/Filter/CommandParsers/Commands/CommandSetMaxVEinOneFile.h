#ifndef SCAFFOLDER_COMMANDSETMAXVEINONEFILE_H
#define SCAFFOLDER_COMMANDSETMAXVEINONEFILE_H

#include "Command.h"

class CommandSetMaxVEinOneFile : public Command {
public:
    void execute(std::string argv, State &state, Filter *filter) override {
        int mv, me;
        std::stringstream ss(argv);
        ss >> mv >> me;
        state.maxVert = mv;
        state.maxEdge = me;
    }
};

#endif //SCAFFOLDER_COMMANDSETMAXVEINONEFILE_H
