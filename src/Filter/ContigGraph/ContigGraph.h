#ifndef SCAFFOLDER_GRAPH_H
#define SCAFFOLDER_GRAPH_H

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <Logger/logger.hpp>
#include <seqan/basic.h>
#include <seqan/gff_io.h>
#include "assert.h"

namespace filter {
    namespace contig_graph {
        class ContigGraph {
        public:
            //Annotation of exon on contig
            struct Exon {
                //Start coordinate of exon on contig
                int b;

                //End coordinate of exon on contig
                int e;

                //Covarage of exon
                double cov;

                //Contig ID
                int id;

                bool operator < (const Exon ex) const {
                    return b < ex.b;
                }
            };

            struct Edge {
                //Edge ID
                int id = 0;

                //start vertex ID
                int from = 0;

                //end vertex ID
                int to = 0;

                //library ID
                int lib = 0;

                //connection weight
                double weight = 0;

                //length of connection(distance between contig), -1 if undefined
                int len = -1;

                //the first coordinate(coordBegin1 <= coordEnd1) of reads alignment on first contig
                int coordBegin1 = 0;

                //the second coordinate(coordBegin1 <= coordEnd1) of reads alignment on first contig
                int coordEnd1 = 0;

                //the first coordinate(coordBegin2 <= coordEnd2) of reads alignment on second contig
                int coordBegin2 = 0;

                //the second coordinate of reads alignment on second contig
                int coordEnd2 = 0;

                // Chromosome name for reference connection
                std::string chr_name = "";

                //extra information
                std::string info = "";

                Edge() {}

                Edge(int id, int from, int to, int lib, int weight, int coordBegin1, int coordEnd1, int coordBegin2,
                     int coordEnd2) :
                        id(id), from(from), to(to), lib(lib), weight(weight), coordBegin1(coordBegin1),
                        coordEnd1(coordEnd1), coordBegin2(coordBegin2), coordEnd2(coordEnd2) {}
            };
            struct Vertex {
                //Vertex ID
                int id;

                //Contig name
                std::string name;

                //Contig length
                int len;

                //IDs of outgoings edges
                std::vector<int> edges;

                //IDs of inconing edges
                std::vector<int> edgesR;

                //List of exons from first strand
                std::vector<Exon> exonsStr1;

                //List of exons from second strand
                std::vector<Exon> exonsStr2;

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

            int mxEdge = 0;

            void addExonBlockFromGffFile(std::string fileName);
        public:
            std::vector<Exon> getExons(int v, int strand);

            //get all edges from vertex v
            std::vector<int> getEdges(int v);

            //get all edges to vertex v
            std::vector<int> getEdgesR(int v);

            //if e: v -> u then to[e] = v
            int getEdgeTo(int e);

            //if e: v -> u then from[e] = u
            int getEdgeFrom(int e);

            //return wieght of edge e
            int getEdgeWeight(int e);

            //return lib for this edge
            int getEdgeLib(int e);

            //get the first coordinate(coordBegin1 <= coordEnd1) of reads alignment on first contig
            int getEdgeCoordB1(int e);

            //get the second coordinate(coordBegin1 <= coordEnd1) of reads alignment on first contig
            int getEdgeCoordE1(int e);

            //get the first coordinate(coordBegin2 <= coordEnd2) of reads alignment on second contig
            int getEdgeCoordB2(int e);

            //get the second coordinate of reads alignment on second contig
            int getEdgeCoordE2(int e);

            //Set the name of chromosome for reference source
            void setEdgeChr(int e, std::string name);

            //return color for this lib
            std::string getLibColor(int l);

            //return name of this lib
            std::string getLibName(int l);

            //return the type of library
            Lib::Type getLibType(int l);

            //get ids of all current sources
            std::vector<int> getLibList();

            //add new vertex with this id, name and len
            int addVertex(int id, std::string name, int len);

            //add to graph new source
            int addLib(std::string color, std::string name, Lib::Type type);

            // get len of contig with id
            int getTargetLen(int id);

            // get name of contig with this id
            std::string getTargetName(int v);

            //get id of contig by contig name
            int getTargetId(std::string name);

            //get count of vertexs
            int getVertexCount();

            //get the count of lib
            int getLibNum();

            //get maximum vertex id
            int getMaxVertId();

            //get extra Info of edge with id e
            std::string getInfo(int e);

            //get Vertex by Id
            Vertex getVertex(int v) {
                targets[v].id = v;
                return targets[v];
            }

            //get list of ids of all available nodes
            std::vector<int> getVertexList();

            //serialize this graph in .gr format in "fileName" file
            void write(std::string fileName);

            //generate ContigGraph from .gr format file
            static ContigGraph read(std::string fileName);

            //add new edge between node with ids v1 and v2 with weight w, length betwwen contigs len and with extra info
            int addEdge(int v, int u, int lib, int w, int b1, int e1, int b2, int e2);

            //change edge weight
            void setWeight(int e, int w);

            //change edge coordinate
            void setCoord(int e, int b1, int e1, int b2, int e2);

            //delte edge with ID e
            void delEdge(int e);

            //delete vertex with ID v
            void delVertex(int v);

            //add annotation from <fileName>
            void addExonBlock(std::string fileName);

            //merge source with id l1 and l2, return new Lib
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
