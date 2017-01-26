#ifndef SCAFFOLDER_COMMANDUPLOADGRAPH_H
#define SCAFFOLDER_COMMANDUPLOADGRAPH_H

#include <Filter/Filters/Filter.h>
#include "Command.h"

//uploadGraph <filename>
class CommandUploadGraph : public Command {
public:
    void execute(std::string argv, State& state, Filter *filter) override;
};


#endif //SCAFFOLDER_COMMANDUPLOADGRAPH_H
