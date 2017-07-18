#ifndef SCAFFOLDER_CLUSTERSTATISTIC_H
#define SCAFFOLDER_CLUSTERSTATISTIC_H

#include "Statistic.h"

namespace filter {
    namespace statistics {
        class ClusterStatistic : Statistic {
        public:
            void calculateStatistic(ContigGraph *filter, std::string coordFile);
        private:
            InfoAboutContigsAlig aligInfo;
            struct clusterInfo {
                int size;
                int minWeight;
                int max2Weight;
                int maxWeight;
                int winWeight;

                bool operator < (clusterInfo cl) {
                    return  (size > cl.size) ||
                            (size == cl.size && maxWeight > cl.maxWeight) ||
                            (size == cl.size && maxWeight == cl.maxWeight && minWeight > cl.minWeight);
                }
            };

            std::vector<clusterInfo> clusters;

            void print();
            void setInfoAboutCluster(filter::ContigGraph *filter, const std::vector<int> &edges,
                    filter::statistics::ClusterStatistic::clusterInfo &clinfo);
        };
    }
}


#endif //SCAFFOLDER_CLUSTERSTATISTIC_H
