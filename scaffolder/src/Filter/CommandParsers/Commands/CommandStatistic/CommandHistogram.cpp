#include <Filter/Statistics/WeightHistogram.h>
#include "CommandHistogram.h"

void CommandHistogram::execute(std::string argv, State &state, Filter *filter) {
    std::stringstream ss(argv);

    int lib, step;
    ss >> lib >> step;

    WeightHistogram::histogram(filter, lib, step);
}
