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
    int libNum;

    vector<vector<int> > graph;
    vector<int> start;
    int startEdgeNum;
    vector<int> to;

    vector<int> edgeLib;
    vector<int> edgeWeight;

    vector<int> targetId;
    vector<string> targetName;
    vector<double> targetCoverage;
    vector<int> targetLen;

    vector<string> libColor;

    vector< vector <int> > edgeIdByVertex;

    string genRandomColor();
public:
    void newLib();
    void filterByEdgeWight(int minEdgeWight);
    void filterByContigLen(int minContigLen);
    void writeGraphDotFormat(string fileName);
    void writeFullGraph();
    void sortEdgeByWight(int v);
    void delEdges(int v, int k);
    vector<int> getEdgesWight(int v);
    void incEdgeWeight(string vName, int u);
    int addVertex(int id, string name, double cov, int len);
    void incVertexCover(int id, double x);
    int getTargetLength(int id);
    int getVertexCount();
    int getLibNum();
};


#endif //SCAFFOLDER_GRAPH_H
