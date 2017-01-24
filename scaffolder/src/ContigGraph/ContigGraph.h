#ifndef SCAFFOLDER_GRAPH_H
#define SCAFFOLDER_GRAPH_H

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include "assert.h"

// Store graph on contigs with several libs
class ContigGraph {
    friend class Serialization;
private:
    std::vector<std::vector<int> > graph; // graph[v][i] = e store edge id from vertex (e: v -> u)
    std::vector<std::vector<int> > graphR; // graph[u][i] = e store edge id to vertex (e: v -> u)
    std::vector<int> to; // if e: v -> u then to[e] = u
    std::vector<int> from; // if e: v -> u then to[e] = v

    std::vector<int> edgeLib; // edgeLib[e] = lib of this edge
    std::vector<int> edgeWeight; // edgeWeight[e] = weight

    std::map<std::string, int> targetId; // return vertex in graph(aka target id) by target name
    std::vector<std::string> targetName; // return target name by target id
    std::vector<int> targetLen; // return target len by target id

    std::vector<std::string> libColor; //return lib color by lib id
    std::vector<std::string> libName; //return lib name by lib id

    std::vector<std::unordered_map<int, int> > vrtsToEdge; // if in last lib e: v -> u then vrtsToEdge[v][u] = e
public:
    void newLib(std::string name, std::string color); //next edge library

    std::vector<int> getEdges(int v); //get all edges from vertex v
    std::vector<int> getEdgesR(int v); //get all edges to vertex v

    int getToVertex(int e); //if e: v -> u then to[e] = v
    int getFromVertex(int e); //if e: v -> u then from[e] = u
    int getEdgeWeight(int e); //return wieght of edge e
    int getEdgeLib(int e); //return lib for this edge
    std::string getLibColor(int l); //return color for this lib
    std::string getLibName(int l); //return name of this lib

    int incEdgeWeight(int v, int u); //increment edge wight between contigs with id v and u

    int addVertex(int id, std::string name, int len); //add new vertex with this id, name and len

    int getTargetLength(int id) const; // get len of contig with id
    std::string getTargetName(int v); // get name of contig with this id
    int getVertexCount(); //get count of vertexs
    int getLibNum(); //get the count of lib

    void write(std::string fileName); //serialize this graph in .gr format in "fileName" file
    static ContigGraph read(std::string fileName); //generate ContigGraph from .gr format file
};
#endif //SCAFFOLDER_GRAPH_H
