#ifndef SCAFFOLDER_FILTERIGNORE_H
#define SCAFFOLDER_FILTERIGNORE_H

#include "Filter.h"

namespace filter {
//filter for ignore some vertexes
    class FilterIgnore : public Filter {
    private:
        std::vector<bool> ignore;
    public:
        FilterIgnore(Filter *filter);

        std::vector<int> getEdges(int v) override;

        std::vector<int> getEdgesR(int v) override;

        std::vector<int> getVertexList() override;

        //handling querys SET_IGNORE with args <first vertex id for ignore> <last vertex id for ignore>
        // will ignore vertex in [first, last)
        //and RESET_IGNORE without args stop ignore all vertexs
        void processQuery(Query query) override;
    };
}


#endif //SCAFFOLDER_FILTERIGNORE_H
