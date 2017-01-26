#ifndef SCAFFOLDER_STATE_H
#define SCAFFOLDER_STATE_H

#include <string>

struct State {
    enum StateName {DEF, LOCAL};
    StateName name = StateName::DEF;
    std::string fileName = "";
    int dist = 0;
};

#endif //SCAFFOLDER_STATE_H
