#ifndef SCAFFOLDER_STATE_H
#define SCAFFOLDER_STATE_H

#include <string>
#include <Filter/Writers/FileValidator/FileValidator.h>

struct State {
    enum StateName {DEF, LOCAL};
    StateName name = StateName::DEF;
    std::string fileName = "";
    int dist = 0;
    FileValidator *validator = new FileValidator;

    ~State() {
        delete validator;
    }
};

#endif //SCAFFOLDER_STATE_H
