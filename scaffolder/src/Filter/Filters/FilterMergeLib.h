#ifndef SCAFFOLDER_FILTERMERGELIB_H
#define SCAFFOLDER_FILTERMERGELIB_H

#include "Filter.h"

namespace filter {
    class FilterMergeLib : public Filter {
    private:
        std::vector<int> newLibNum;
        std::vector<std::string> newLibName;
        std::vector<double> wc;
    public:
        FilterMergeLib(Filter *filter);

        std::vector<int> getLibList() override;

        std::vector<int> getEdges(int v) override;

        std::vector<int> getEdgesR(int v) override;

        int getEdgeWeight(int e) override;

        int getEdgeLib(int e) override;

        std::string getLibName(int l) override;

        std::string getLibColor(int l) override;

        void processQuery(Query query) override;

        void update(Filter *pFilter);
    };
}

#endif //SCAFFOLDER_FILTERMERGELIB_H
