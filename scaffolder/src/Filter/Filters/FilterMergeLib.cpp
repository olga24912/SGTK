#include <algorithm>
#include <iostream>
#include "FilterMergeLib.h"

namespace filter {
    std::vector<int> FilterMergeLib::getLibList() {
        TRACE("getLibList");
        std::vector<int> libList;
        for (int i = 0; i < (int) newLibNum.size(); ++i) {
            libList.push_back(newLibNum[i]);
        }

        std::sort(libList.begin(), libList.end());

        libList.resize(std::unique(libList.begin(), libList.end()) - libList.begin());

        return libList;
    }

    std::vector<int> FilterMergeLib::getEdges(int v) {
        TRACE("getEdges v=" << v);
        std::vector<int> edges = Filter::getEdges(v);
        std::vector<std::pair<std::pair<int, int>, int>> edgesInfo;
        for (int i : edges) {
            edgesInfo.push_back(std::make_pair(std::make_pair(getEdgeTo(i),
                                                              FilterMergeLib::getEdgeLib(i)), i));
        }

        sort(edgesInfo.begin(), edgesInfo.end());

        std::vector<int> res;
        for (int i = 0; i < (int) edgesInfo.size(); ++i) {
            if (i == 0 || edgesInfo[i - 1].first != edgesInfo[i].first) {
                res.push_back(edgesInfo[i].second);
            }
        }

        return res;
    }

    std::vector<int> FilterMergeLib::getEdgesR(int v) {
        TRACE("getEdgesR v=" << v);
        std::vector<int> edges = Filter::getEdgesR(v);
        std::vector<std::pair<std::pair<int, int>, int>> edgesInfo;
        for (int i : edges) {
            edgesInfo.push_back(std::make_pair(std::make_pair(getEdgeFrom(i),
                                                              FilterMergeLib::getEdgeLib(i)), i));
        }

        sort(edgesInfo.begin(), edgesInfo.end());

        std::vector<int> res;
        for (int i = 0; i < (int) edgesInfo.size(); ++i) {
            if (i == 0 || edgesInfo[i - 1].first != edgesInfo[i].first) {
                res.push_back(edgesInfo[i].second);
            }
        }

        return res;
    }

    int FilterMergeLib::getEdgeWeight(int e) {
        TRACE("getEdgeWeight e=" << e);
        int v = Filter::getEdgeFrom(e), u = Filter::getEdgeTo(e);

        int w = 0;
        std::vector<int> edges = Filter::getEdges(v);
        for (int i : edges) {
            if (getEdgeLib(i) == getEdgeLib(e) && getEdgeTo(i) == u) {
                w += Filter::getEdgeWeight(i) * wc[Filter::getEdgeLib(i)];
            }
        }

        return w;
    }

    int FilterMergeLib::getEdgeLib(int e) {
        TRACE("getEdgeLib e=" << e);
        return newLibNum[Filter::getEdgeLib(e)];
    }

    std::string FilterMergeLib::getLibName(int l) {
        TRACE("getLibName l=" << l);
        return newLibName[newLibNum[l]];
    }

    std::string FilterMergeLib::getLibColor(int l) {
        TRACE("getLibColor l=" << l);
        return Filter::getLibColor(newLibNum[l]);
    }

    void FilterMergeLib::processQuery(Query query) {
        TRACE("process query");
        if (query.type != query.MERGE_LIB) {
            Filter::processQuery(query);
            update(subfilter);
        } else {
            std::stringstream ss(query.argv);
            int n1, n2;
            std::string name;
            ss >> n1 >> n2 >> name;

            double w1 = 1, w2 = 1;
            ss >> w1 >> w2;

            TRACE("query merge_lib n1=" << n1 << " n2=" << n2 << " name=" << name << " with weight = " << w1 << " " << w2);


            for (int i = 0; i < newLibNum.size(); ++i) {
                if (newLibNum[i] == n1) {
                    wc[i] *= w1;
                }
                if (newLibNum[i] == n2) {
                    wc[i] *= w2;
                }
            }

            int k = n1;
            n1 = std::min(n1, n2);
            n2 = std::max(k, n2);

            newLibName[n1] = name;

            for (int i = 0; i < newLibNum.size(); ++i) {
                if (newLibNum[i] == n2) {
                    newLibNum[i] = n1;
                }
            }
        }
    }

    FilterMergeLib::FilterMergeLib(Filter *filter) : Filter(filter) {
        update(filter);
    }

    void FilterMergeLib::update(Filter *filter) {
        TRACE("update");
        int libCnt = (int) filter->getLibList().size();
        wc.resize(0);
        wc.resize(libCnt, 1);
        newLibNum.resize(libCnt);
        newLibName.resize(libCnt);
        for (int i = 0; i < libCnt; ++i) {
            newLibNum[i] = i;
            newLibName[i] = filter->getLibName(i);
        }
    }
}
