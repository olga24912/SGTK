#include <Filter/Statistics/cntRightConnectionStatistic.h>
#include "CommandCorrectConnectionStatistic.h"

void CommandCorrectConnectionStatistic::execute(std::string argv, State &state, Filter *filter) {
    std::stringstream ss(argv);
    std::string fileName;
    int libNum;
    ss >> fileName >> libNum;

    cntRightConnectionStatistic cs;
    cs.calculateStatistic(filter, fileName, libNum);
}
