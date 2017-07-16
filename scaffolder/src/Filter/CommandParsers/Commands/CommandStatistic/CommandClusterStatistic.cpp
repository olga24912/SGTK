#include <Filter/Statistics/ClusterStatistic.h>
#include "CommandClusterStatistic.h"

void filter::commands::CommandClusterStatistic::execute(std::string argv, filter::commands::State &state,
                                                        filter::Filter *filter) {
    std::stringstream ss(argv);
    std::string fileName;
    ss >> fileName;

    INFO("Cluster statistic fileName=" << fileName);

    statistics::ClusterStatistic cs;
    cs.calculateStatistic(filter, fileName);
}
