#ifndef SCAFFOLDER_COMMANDSETFV_H
#define SCAFFOLDER_COMMANDSETFV_H

#include <Filter/CommandParsers/Commands/Command.h>

class CommandSetFV : public Command {
public:
    void execute(std::string argv, State &state, Filter *filter) override final {
        delete state.validator;
        setFV(state, argv);
    }

protected:
    virtual void setFV(State &state, std::string argv) {
        state.validator = new FileValidator;
    }
};

#endif //SCAFFOLDER_COMMANDSETFV_H
