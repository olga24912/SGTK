#include "ClusterStatistic.h"

namespace filter {
    namespace statistics {
        void filter::statistics::ClusterStatistic::calculateStatistic(ContigGraph *graph, std::string coordFile) {
            aligInfo.parseCoordFile(graph, coordFile);

            std::vector<int> vert = graph->getVertexList();

            for (int v : vert) {
                std::vector<int> edges = graph->getEdges(v);

                clusterInfo clinfo;

                setInfoAboutCluster(graph, edges, clinfo);
            }

            std::sort(clusters.begin(), clusters.end());

            print();
        }

        void ClusterStatistic::setInfoAboutCluster(ContigGraph *filter, const std::vector<int> &edges,
                                                   filter::statistics::ClusterStatistic::clusterInfo &clinfo) {
            clinfo.size = edges.size();
            clinfo.winWeight = 0;
            clinfo.minWeight = 1e9;
            clinfo.maxWeight = 0;
            std::vector<int> weights;
            for (int e : edges) {
                int ww = filter->getEdgeWeight(e);
                weights.push_back(ww);
                if (aligInfo.isCorrectEdge(filter, e) == filter::alig_info::InfoAboutContigsAlig::OK) {
                    clinfo.winWeight = std::max(clinfo.winWeight, ww);
                }
            }

            if (clinfo.size == 0) return;

            std::sort(weights.begin(), weights.end());
            clinfo.minWeight = weights[0];
            clinfo.maxWeight = weights[weights.size() - 1];
            clinfo.max2Weight = clinfo.maxWeight;
            if (clinfo.size > 1) {
                clinfo.max2Weight = weights[weights.size() - 2];
            }

            if (clinfo.size > 0) {
                this->clusters.push_back(clinfo);
            }
        }

        void ClusterStatistic::print() {
            for (int i = 0; i < (int)clusters.size(); ++i) {
                std::cout << clusters[i].size << "\t"  << clusters[i].minWeight << "\t"
                          << clusters[i].max2Weight << "\t" << clusters[i].maxWeight
                          << "\t" << clusters[i].winWeight << "\t"  << (clusters[i].winWeight != 0) << "\n";
            }
        }
    }
}