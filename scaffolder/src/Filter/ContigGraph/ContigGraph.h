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

namespace filter {
    namespace contig_graph {
        // Store graph on contigs with several libs
        class ContigGraph {
        public:
            struct Edge {
                int id = 0;
                int from = 0;
                int to = 0;
                int lib = 0;
                int weight = 0;
                int coordBegin1 = 0;
                int coordEnd1 = 0;
                int coordBegin2 = 0;
                int coordEnd2 = 0;
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
                std::vector<int> edges;
                std::vector<int> edgesR;

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

        private:
            std::unordered_map<std::string, int> targetId; // return vertex in graph(aka target id) by target name

            std::map<int, Vertex> targets;
            std::map<int, Lib> libs;
            std::map<int, Edge> edges;

            int edgeCnt;
            int vertCnt;
            int libCnt;
        public:
            std::vector<int> getEdges(int v); //get all edges from vertex v
            std::vector<int> getEdgesR(int v); //get all edges to vertex v

            int getEdgeTo(int e); //if e: v -> u then to[e] = v
            int getEdgeFrom(int e); //if e: v -> u then from[e] = u
            int getEdgeWeight(int e); //return wieght of edge e
            int getEdgeLib(int e); //return lib for this edge
            int getEdgeCoordB1(int e);
            int getEdgeCoordE1(int e);
            int getEdgeCoordB2(int e);
            int getEdgeCoordE2(int e);

            void setEdgeChr(int e, std::string name);

            std::string getLibColor(int l); //return color for this lib
            std::string getLibName(int l); //return name of this lib
            Lib::Type getLibType(int l);

            std::vector<int> getLibList();

            int addVertex(int id, std::string name, int len); //add new vertex with this id, name and len

            int getTargetLen(int id); // get len of contig with id
            std::string getTargetName(int v); // get name of contig with this id
            int getTargetId(std::string name); //get id of contig by contig name
            int getVertexCount(); //get count of vertexs
            int getLibNum(); //get the count of lib
            int getMaxVertId();
            std::string getInfo(int e);

            std::vector<int> getVertexList();

            void write(std::string fileName); //serialize this graph in .gr format in "fileName" file
            static ContigGraph read(std::string fileName); //generate ContigGraph from .gr format file

            int addEdge(int v, int u, int lib, int w, int b1, int e1, int b2, int e2);
            void setWeight(int e, int w);
            void delEdge(int e);
            void delVertex(int v);

            Lib mergeLib(int l1, int l2, std::string lib_name, double w1, double w2);
        private:
            DECL_LOGGER("ContigGraph");

            const static int COORD_DIST = 50;

            void mergeEdges(Vertex v, int l1);

            void mergeEdge(int &e1, int &e2);
            void delLib(int l);

            bool intersec(int b1, int e1, int b2, int e2);
        };
    }
}
#endif //SCAFFOLDER_GRAPH_H
