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
        class ContigGraph {
        public:
            const int MAX_EDGES_CNT=5*1e5;
            struct Edge {
                //Edge ID
                int id;

                //start vertex ID
                int from;

                //end vertex ID
                int to;

                //library ID
                int lib;

                //connection weight
                double weight;

                //length of connection(distance between contig), -1 if undefined
                int len = -1;

                //the first coordinate(coordBegin1 <= coordEnd1) of reads alignment on first contig
                int coordBegin1;

                //the second coordinate(coordBegin1 <= coordEnd1) of reads alignment on first contig
                int coordEnd1;

                //the first coordinate(coordBegin2 <= coordEnd2) of reads alignment on second contig
                int coordBegin2;

                //the second coordinate of reads alignment on second contig
                int coordEnd2;

                //extra information
                std::string info = "";

                Edge() {}

                Edge(int id, int from, int to, int lib, int weight, int coordBegin1, int coordEnd1, int coordBegin2,
                     int coordEnd2) :
                        id(id), from(from), to(to), lib(lib), weight(weight), coordBegin1(coordBegin1),
                        coordEnd1(coordEnd1), coordBegin2(coordBegin2), coordEnd2(coordEnd2) {}


                Edge(int id, int from, int to, int lib, double weight, int coordBegin1, int coordEnd1, int coordBegin2,
                     int coordEnd2) :
                        id(id), from(from), to(to), lib(lib), weight(weight), coordBegin1(coordBegin1),
                        coordEnd1(coordEnd1), coordBegin2(coordBegin2), coordEnd2(coordEnd2) {}
            };

            struct Vertex {
                //vertex ID
                int id;

                //name of contig
                std::string name;

                //contig length
                int len;

                Vertex() {}

                Vertex(int id, std::string name, int len) : id(id), name(name), len(len) {}
            };

            struct Lib {
                static const int typeCnt = 12;
                enum Type {
                    REF, DNA_PAIR, RNA_PAIR, RNA_SPLIT_50, RNA_SPLIT_30, SCAFF, CONNECTION, MATE_PAIR, LONG, FASTG, GFA, GFA2
                };
                static const std::string typeToStr[];

                //TODO: delete
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
            static const int maxClusterSize = 1000;

            // graph[v][i] = e store edge id from vertex (e: v -> u)
            std::vector<std::vector<int> > graph;

            // graph[u][i] = e store edge id to vertex (e: v -> u)
            std::vector<std::vector<int> > graphR;

            // return vertex in graph(aka target id) by target name
            std::unordered_map<std::string, int> targetId;

            std::vector<Vertex> targets;
            std::vector<Lib> libs;
            std::vector<Edge> edges;

            void filterGraph();
        public:
            //next edge library
            void newLib(std::string name, std::string color, Lib::Type type);

            //get all edges ids from vertex v
            std::vector<int> getEdges(int v);

            //get all edges ids to vertex v
            std::vector<int> getEdgesR(int v);

            //get all edges between vertex v and vertex u
            std::vector<Edge> getEdgesBetween(int v, int u);

            //return color for this library
            std::string getLibColor(int l);

            //return type of library
            Lib::Type getLibType(int l);

            //increment edge wight between contigs with id v and u
            int incEdgeWeight(int v, int u, int cb1, int ce1, int cb2, int ce2);

            //add new vertex with this id, name and len
            int addVertex(int id, std::string name, int len);

            // get len of contig with id
            int getTargetLen(int id) const;

            //get count of vertexs
            int getVertexCount();

            //get the count of lib
            int getLibNum();

            //get vertex ID by contig name
            int getTargetId(std::string name);

            //get contig name by vertex id
            std::string getTargetName(int id) const;

            //add new edge between nodes with ids v1 and v2 with coordinates of reads alignment c1 and c2
            int addEdge(int v1, int v2, std::pair<int, int> c1, std::pair<int, int> c2);

            //add new edge between node with ids v1 and v2 with weight w, length betwwen contigs len and with extra info
            int addEdge (int v1, int v2, double w, int len=-1, std::string info = "");

            //increment edge weight with ID e, and update coordinate
            void incEdge(int e, std::pair<int, int> c1, std::pair<int, int> c2);

            //serialize this graph in .gr format in "fileName" file
            void write(std::string fileName);

            //generate ContigGraph from .gr format file
            static ContigGraph read(std::string fileName);
        private:
            DECL_LOGGER("ContigGraph");
        };
    }
}
#endif //SCAFFOLDER_GRAPH_H
