#ifndef SCAFFOLDER_CLUSTERSTATISTIC_H
#define SCAFFOLDER_CLUSTERSTATISTIC_H

#include <ContigGraph/ContigGraph.h>
#include "Statistic.h"

namespace filter {
    namespace statistics {
        using namespace contig_graph;
        class ClusterStatistic : Statistic {
        public:
            void calculateStatistic(ContigGraph *graph, std::string coordFile);
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
            void setInfoAboutCluster(ContigGraph *filter, const std::vector<int> &edges,
                    filter::statistics::ClusterStatistic::clusterInfo &clinfo);
        };
    }
}


#endif //SCAFFOLDER_CLUSTERSTATISTIC_H
