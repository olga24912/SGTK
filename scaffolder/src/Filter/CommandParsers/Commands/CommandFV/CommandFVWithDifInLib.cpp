#include <Filter/Writers/FileValidator/ValidatorWithDifInLib.h>
#include "CommandFVWithDifInLib.h"

void CommandFVWithDifInLib::setFV(State &state, std::string argv) {
    int libNum;
    std::vector<int> libs;
    std::stringstream ss(argv);
    while (ss >> libNum) {
        libs.push_back(libNum);
    }

    state.validator = new ValidatorWithDifInLib(libs);
}
