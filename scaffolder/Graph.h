//
// Created by olga on 09.10.16.
//

#ifndef SCAFFOLDER_GRAPH_H
#define SCAFFOLDER_GRAPH_H

#include <bits/stdc++.h>
#include <seqan/bam_io.h>
#include <seqan/graph_types.h>

using namespace std;
using namespace seqan;

class Graph {
private:
    int libNum = -1;

    int minContigLen;

    vector<vector<int> > graph;
    vector<int> start;
    int startEdgeNum = 0;
    vector<int> to;

    vector<int> edgeLib;
    vector<int> edgeWeight;

    map<string, int> targetId;
    vector<int> vById;
    vector<int> idByV;
    vector<string> targetName;
    vector<double> targetCoverage;
    vector<int> targetLen;

    vector<string> libColor;

    vector< vector <int> > edgeIdByVertexes;

    string genRandomColor();
public:
    void newLib();
    void filterByEdgeWeight(int minEdgeWeight);
    void filterByContigLen(int minContigLen);
    void writeGraphDotFormat(string fileName);
    void sortEdgeByWeight(int v);
    void delEdges(int v, int k);
    vector<int> getEdgesWeight(int v);
    void incEdgeWeight(int vId, int uId);
    int addVertex(int id, string name, double cov, int len);
    void incTargetCover(int id, double x);
    int getTargetLength(int id);
    int getVertexCount();
    int getLibNum();
};


#endif //SCAFFOLDER_GRAPH_H
