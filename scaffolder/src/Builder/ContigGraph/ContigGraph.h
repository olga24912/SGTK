#ifndef SCAFFOLDER_GRAPH_H
#define SCAFFOLDER_GRAPH_H

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <Logger/logger.hpp>
#include "assert.h"

namespace builder {
    namespace contig_graph {
// Store graph on contigs with several libs
        class ContigGraph {
        public:
            struct Edge {
                int id;
                int from;
                int to;
                int lib;
                int weight;
                int coordBegin1;
                int coordEnd1;
                int coordBegin2;
                int coordEnd2;
                std::string chr_name = "";

                Edge() {}

                Edge(int id, int from, int to, int lib, int weight, int coordBegin1, int coordEnd1, int coordBegin2,
                     int coordEnd2) :
                        id(id), from(from), to(to), lib(lib), weight(weight), coordBegin1(coordBegin1),
                        coordEnd1(coordEnd1), coordBegin2(coordBegin2), coordEnd2(coordEnd2) {}
            };

            struct Vertex {
                int id;
                std::string name;
                int len;

                Vertex() {}

                Vertex(int id, std::string name, int len) : id(id), name(name), len(len) {}
            };

            struct Lib {
                static const int typeCnt = 6;
                enum Type {
                    REF, DNA_PAIR, RNA_PAIR, RNA_SPLIT_50, RNA_SPLIT_30, SCAFF
                };
                static const std::string typeToStr[];
                std::string color;
                std::string name;
                Type type;

                Lib() {}

                Lib(std::string color, std::string name, Type type) : color(color), name(name), type(type) {}

                Lib(std::string color, std::string name, std::string type) : color(color), name(name) {
                    for (int i = 0; i < typeCnt; ++i) {
                        if (typeToStr[i] == type) {
                            this->type = (Type) i;
                        }
                    }
                }
            };

            int getTargetId(std::string name);

        private:
            static const int maxClusterSize = 1000;
            std::vector<std::vector<int> > graph; // graph[v][i] = e store edge id from vertex (e: v -> u)
            std::vector<std::vector<int> > graphR; // graph[u][i] = e store edge id to vertex (e: v -> u)
            std::map<std::string, int> targetId; // return vertex in graph(aka target id) by target name

            std::vector<Vertex> targets;
            std::vector<Lib> libs;
            std::vector<Edge> edges;

            std::vector<std::unordered_map<int, std::vector<int> > > vrtsToEdge; // if in last lib e: v -> u then vrtsToEdge[v][u] = e
        public:
            void newLib(std::string name, std::string color, Lib::Type type); //next edge library

            std::vector<int> getEdges(int v); //get all edges from vertex v
            std::vector<int> getEdgesR(int v); //get all edges to vertex v
            int getEdgeTo(int e); //if e: v -> u then to[e] = v
            int getEdgeFrom(int e); //if e: v -> u then from[e] = u
            int getEdgeWeight(int e); //return wieght of edge e
            int getEdgeLib(int e); //return lib for this edge
            void setEdgeChr(int e, std::string name);

            std::string getLibColor(int l); //return color for this lib
            Lib::Type getLibType(int l);
            void setLib(Lib::Type type);

            int incEdgeWeight(int v, int u, int cb1, int ce1, int cb2,
                              int ce2); //increment edge wight between contigs with id v and u
            int addVertex(int id, std::string name, int len); //add new vertex with this id, name and len
            int getTargetLen(int id) const; // get len of contig with id
            int getVertexCount(); //get count of vertexs
            int getLibNum(); //get the count of lib

            void write(std::string fileName); //serialize this graph in .gr format in "fileName" file
            static ContigGraph read(std::string fileName); //generate ContigGraph from .gr format file
        private:
            DECL_LOGGER("ContigGraph");
        };
    }
}
#endif //SCAFFOLDER_GRAPH_H
