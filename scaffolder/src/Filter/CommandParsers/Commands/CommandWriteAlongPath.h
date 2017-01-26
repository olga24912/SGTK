#ifndef SCAFFOLDER_COMMANDWRITEALONGPATH_H
#define SCAFFOLDER_COMMANDWRITEALONGPATH_H

#include "CommandWrite.h"
#include <Filter/Writers/WriteAlongPath.h>

class CommandWriteAlongPath : public CommandWrite {
public:
    void writeGraph(std::string argv, State &state, Filter *filter) override;
};


#endif //SCAFFOLDER_COMMANDWRITEALONGPATH_H
