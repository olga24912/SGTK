#ifndef SCAFFOLDER_COMMANDWRITEFULL_H
#define SCAFFOLDER_COMMANDWRITEFULL_H

#include "CommandWrite.h"
#include <Filter/Writers/WriteFullGraph.h>

class CommandWriteFull : public CommandWrite {
public:
    void writeGraph(std::string argv, State& state, Filter *filter) override;
};

#endif //SCAFFOLDER_COMMANDWRITEFULL_H
