#include <Filter/Statistics/DifStatistic.h>
#include "CommandDifStatistic.h"

void CommandDifStatistic::execute(std::string argv, State &state, Filter *filter) {
    std::stringstream ss(argv);
    std::string fileName;
    int libNum;
    int step;
    int mxVal;
    ss >> fileName >> libNum >> step >> mxVal;

    INFO("Dif statistic fileName=" << fileName << " libNum=" << libNum << " step=" << step << " mxVal=" << mxVal);

    DifStatistic cs;
    cs.calculateStatistic(filter, fileName, libNum, step, mxVal);
}
