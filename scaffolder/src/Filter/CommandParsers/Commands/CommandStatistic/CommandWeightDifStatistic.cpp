#include <Filter/Statistics/WeightDifStatistic.h>
#include "CommandWeightDifStatistic.h"

void CommandWeightDifStatistic::execute(std::string argv, State &state, Filter *filter) {
    std::stringstream ss(argv);
    std::string fileName;
    int libNum;
    int stepW, stepD;
    int mxW, mxD;
    ss >> fileName >> libNum >> stepW >> mxW >> stepD >> mxD;

    WeightDifStatistic cs;
    cs.calculateStatistic(filter, fileName, libNum, stepW, mxW, stepD, mxD);
}
