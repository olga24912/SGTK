//
// Created by olga on 09.10.16.
//

#ifndef SCAFFOLDER_GRAPH_H
#define SCAFFOLDER_GRAPH_H

#include <bits/stdc++.h>
#include <seqan/bam_io.h>
#include <seqan/graph_types.h>
#include "GraphUtils.h"

using namespace std;
using namespace seqan;

// Store graph on contigs with several libs
class ContigGraph {
    friend class Serialization;
private:
    vector<vector<int> > graph; // graph[v][i] = e store edge id from vertex (e: v -> u)
    vector<vector<int> > graphR; // graph[u][i] = e store edge id to vertex (e: v -> u)
    vector<int> to; // if e: v -> u then to[e] = u
    vector<int> from; // if e: v -> u then to[e] = v

    vector<int> edgeLib; // edgeLib[e] = lib of this edge
    vector<int> edgeWeight; // edgeWeight[e] = weight

    map<string, int> targetId; // return vertex in graph(aka target id) by target name
    vector<string> targetName; // return target name by target id
    vector<int> targetLen; // return target len by target id

    vector<string> libColor; //return lib color by lib id
    vector<string> libName; //return lib name by lib id

    vector<unordered_map<int, int> > vrtsToEdge; // if in last lib e: v -> u then vrtsToEdge[v][u] = e
public:
    void newLib(string name, string color); //next edge library

    vector<int> getEdges(int v); //get all edges from vertex v
    vector<int> getEdgesR(int v); //get all edges to vertex v

    int getToVertex(int e); //if e: v -> u then to[e] = v
    int getFromVertex(int e); //if e: v -> u then from[e] = u

    int incEdgeWeight(int v, int u); //increment edge wight between contigs with id v and u

    int addVertex(int id, string name, int len); //add new vertex with this id, name and len

    int getTargetLength(int id) const; // get len of contig with id.
    int getVertexCount(); //get count of vertexs
    int getLibNum(); //get the count of lib
};
#endif //SCAFFOLDER_GRAPH_H
