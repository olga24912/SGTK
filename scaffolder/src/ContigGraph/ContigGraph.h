//
// Created by olga on 09.10.16.
//

#ifndef SCAFFOLDER_GRAPH_H
#define SCAFFOLDER_GRAPH_H

#include <bits/stdc++.h>
#include <seqan/bam_io.h>
#include <seqan/graph_types.h>
#include "Tools/Utils.h"

using namespace std;
using namespace seqan;

/*
 * class with contis graph
 */

class ContigGraph {
    friend class ContigGraphPrinter;
    friend class Serialization;
private:
    int minContigLen;

    vector<vector<int> > graph;
    vector<vector<int> > graphR;
    vector<int> start;
    int startEdgeNum = 0;
    vector<int> to;
    vector<int> from;

    vector<int> edgeLib;
    vector<int> edgeWeight;

    map<string, int> targetId;
    vector<int> vById;
    vector<int> idByV;
    vector<string> targetName;
    vector<double> targetCoverage;
    vector<int> targetLen;

    vector<string> libColor;
    vector<string> libName;
    vector<int> libMinEdgeWight;

    vector< vector <int> > edgeIdByVertexes;

    vector<bool> ignore;

    string genRandomColor();
public:
    void newLib(); //next edge library

    void filterByEdgeWeight(int minEdgeWeight); //delete edge with small wight

    void filterByContigLen(int minContigLen);//fileter vertex with small len

    void sortEdgeByWeight(int v); // sort edge by weight for this vertex

    void delEdges(int v, int k); // delete last k edges for vertex v

    vector<int> getEdgesWeight(int v); //get out edges weight for vertex v.

    vector<int> getEdges(int v);

    vector<int> getEdgesR(int v);

    int getToVertex(int e);

    int getFromVertex(int e);

    int incEdgeWeight(int vId, int uId); //increment edge wight between contigs with id vId and uId

    int addVertex(int id, string name, double cov, int len); //add new vertex with this id, name, coverage and len

    void incTargetCover(int id, double x); //increment coverage for contig with this id on x.

    int getTargetLength(int id)const; // get len of contig with id.

    int getVertexCount(); //get count of vertexs

    int getLibNum(); //get the count of lib

    void setLibName(string s);

    void setMinEdgeWeightForLib(int libNum, int minWeight);

    bool isGoodEdge(int e);

    bool isGoodVertex(int v);

    void setColor(string color);

    void setIgnore(int v) {
        ignore[v] = ignore[v] ^ 1;
    }

    bool isIgnore(int v) {
        return ignore[v];
    }
};


#endif //SCAFFOLDER_GRAPH_H
