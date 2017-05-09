#include <Filter/Statistics/TwoLibStatistic.h>
#include "CommandTwoLibStatistic.h"

void CommandTwoLibStatistic::execute(std::string argv, State &state, Filter *filter) {
    std::stringstream ss(argv);
    std::string fileName;
    int libNum1;
    int step1;
    int mxVal1;
    int libNum2;
    int step2;
    int mxVal2;
    ss >> fileName >> libNum1 >> step1 >> mxVal1 >> libNum2 >> step2 >> mxVal2;

    TwoLibStatistic cs;
    cs.calculateStatistic(filter, fileName, libNum1, step1, mxVal1, libNum2, step2, mxVal2);
}
