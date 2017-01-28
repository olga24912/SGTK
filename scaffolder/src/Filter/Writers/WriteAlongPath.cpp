//
// Created by olga on 26.01.17.
//

#include "WriteAlongPath.h"

WriteAlongPath::WriteAlongPath(std::string fileName, int libId, int dist, int minSize, Filter *filter1) :
        fileName(fileName), libId(libId), dist(dist), minSize(minSize), Writer(filter1) {}

void WriteAlongPath::write() {
    int n = (filter->getVertexCount());
    int *col = new int[n];
    int cur = searcher.findComponent(col);

    std::vector< std::vector<int> > parts(cur);
    for (int i = 0; i < n; ++i) {
        if (col[i] == 0) continue;
        parts[col[i]].push_back(i);
    }

    for (int i = 1; i < cur; ++i) {
        if (parts[i].size() < minSize) continue;
        std::string fn = fileName;
        std::stringstream ss;
        ss << i;
        fn += ss.str();

        std::vector<int> res;
        for (int j = 0; j < (int)parts[i].size(); ++j) {
            std::vector<int> local = searcher.findVertInLocalArea(parts[i][j], dist);
            for (int g = 0; g < (int)local.size(); ++g) {
                res.push_back(local[g]);
            }
        }

        sort(res.begin(), res.end());
        res.resize(unique(res.begin(), res.end()) - res.begin());

        dotWriter.writeVertexSet(res, fn);
    }
    delete col;
}