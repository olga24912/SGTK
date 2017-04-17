#ifndef SCAFFOLDER_COMMANDWRITEBIGCOMP_H
#define SCAFFOLDER_COMMANDWRITEBIGCOMP_H

#include "CommandWrite.h"
#include <Filter/Writers/WriteBigComponent.h>

class CommandWriteBigComp : public CommandWrite {
public:
    void writeGraph(std::string argv, State &state, Filter *filter) override;
};


#endif //SCAFFOLDER_COMMANDWRITEBIGCOMP_H
