#include <Filter/Statistics/WeightStatistic.h>
#include "CommandWeightStatistic.h"

void CommandWeightStatistic::execute(std::string argv, State &state, Filter *filter) {
    std::stringstream ss(argv);
    std::string fileName;
    int libNum;
    int step;
    int mxVal;
    ss >> fileName >> libNum >> step >> mxVal;

    WeightStatistic cs;
    cs.calculateStatistic(filter, fileName, libNum, step, mxVal);
}
