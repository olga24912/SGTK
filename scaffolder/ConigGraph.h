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

/*
 * class with contis graph
 */
class ConigGraph {
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
    /*
     * next edge library
     */
    void newLib();

    /*
     * delete edge with small wight
     */
    void filterByEdgeWeight(int minEdgeWeight);

    /*
     * fileter vertex with small len
     */
    void filterByContigLen(int minContigLen);

    /*
     * write info about conig graph in dot format
     */
    void writeGraphDotFormat(string fileName);

    /*
     * sort edge by weight for this vertex
     */
    void sortEdgeByWeight(int v);

    /*
     * delete last k edges for vertex v
     */
    void delEdges(int v, int k);

    /*
     * get out edges weight for vertex v.
     */
    vector<int> getEdgesWeight(int v);

    /*
     * increment edge wight between contigs with id vId and uId
     */
    void incEdgeWeight(int vId, int uId);

    /*
     * add new vertex with this id, name, coverage and len
     */
    int addVertex(int id, string name, double cov, int len);

    /*
     * increment coverage for contig with this id on x.
     */
    void incTargetCover(int id, double x);

    /*
     * get len of contig with id.
     */
    int getTargetLength(int id)const;

    /*
     * get count of vertexs
     */
    int getVertexCount();

    /*
     * get the count of lib
     */
    int getLibNum();
};


#endif //SCAFFOLDER_GRAPH_H
