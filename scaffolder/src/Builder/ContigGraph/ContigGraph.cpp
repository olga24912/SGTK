#include <sstream>
#include "ContigGraph.h"
#include "cmath"

namespace builder {
    namespace contig_graph {
        const std::string ContigGraph::Lib::typeToStr[] = {"REF", "DNA_PAIR", "RNA_PAIR", "RNA_SPLIT_50",
                                                           "RNA_SPLIT_30",
                                                           "SCAFF", "CONNECTION"};

        int ContigGraph::getLibNum() {
            TRACE("get lib num: " << (int) libs.size());
            return (int) libs.size();
        }

        int ContigGraph::getVertexCount() {
            TRACE("get vertex count: " << (int) graph.size());
            return (int) graph.size();
        }

        int ContigGraph::getTargetLen(int v) const {
            TRACE("get target " << v << " length: " << targets[v].len);
            return targets[v].len;
        }

        int ContigGraph::addVertex(int id, std::string name, int len) {
            assert(id == (int) graph.size());
            DEBUG("add vertex id = " << id << "  vertex name = " << name << "  vertex len = " << len);

            graph.push_back(std::vector<int>());
            graphR.push_back(std::vector<int>());
            targetId[name] = id;

            targets.resize((size_t) id + 1);
            targets[id].name = name;
            targets[id].len = len;

            return id;
        }

        int ContigGraph::incEdgeWeight(int v, int u, int cb1, int ce1, int cb2, int ce2) {
            assert(libs.size() > 0);

            DEBUG("increment edge Weight between v=" << v << "(" << cb1 << "-" << ce1 << ") u=" << u << "(" << cb2
                                                     << "-"
                                                     << ce2 << ")");

            int e = -1;
            std::vector<int> candidates = graph[v];

            for (int ec : candidates) {
                if (std::fabs(cb1 - edges[ec].coordBegin1) < maxClusterSize &&
                    std::fabs(cb2 - edges[ec].coordBegin2) < maxClusterSize) {
                    e = ec;
                }
            }

            if (e == -1) {
                e = (int) edges.size();
                edges.push_back(Edge(e, v, u, (int) libs.size() - 1, 0, cb1, ce1, cb2, ce2));
                graph[v].push_back(e);
                graphR[u].push_back(e);
            }

            edges[e].weight += 1;
            edges[e].coordBegin1 = std::min(edges[e].coordBegin1, cb1);
            edges[e].coordEnd1 = std::max(edges[e].coordEnd1, ce1);
            edges[e].coordBegin2 = std::min(edges[e].coordBegin2, cb2);
            edges[e].coordEnd2 = std::max(edges[e].coordEnd2, ce2);

            DEBUG("inc edge id=" << e << "coords =(" << edges[e].coordBegin1 << "-" << edges[e].coordEnd1
                                 << "; " << edges[e].coordBegin2 << "-" << edges[e].coordEnd2 << ")");
            return e;
        }

        void ContigGraph::newLib(std::string name, std::string color, Lib::Type type) {
            DEBUG("new lib with name=" << name << " color=" << color << " type=" << Lib::typeToStr[type]);
            libs.push_back(Lib(color, name, type));
        }

        std::vector<int> ContigGraph::getEdges(int v) {
            TRACE("getEdges v=" << v);
            return graph[v];
        }

        std::vector<int> ContigGraph::getEdgesR(int v) {
            TRACE("getEdgesR v=" << v);
            return graphR[v];
        }

        void ContigGraph::write(std::string fileName) {
            INFO("start write graph to " << fileName);
            std::ofstream out(fileName);
            DEBUG("libs num=" << libs.size() << " vertex num=" << graph.size() << " edegs num=" << edges.size());

            out << libs.size() << "\n";
            for (int i = 0; i < (int) libs.size(); ++i) {
                out << "l " << i << " " << libs[i].color << " " << libs[i].name << " " << Lib::typeToStr[libs[i].type]
                    << "\n";
            }
            out << graph.size() << "\n";
            for (int i = 0; i < (int) graph.size(); ++i) {
                out << "v " << i << " " << targets[i].name << " " << targets[i].len << "\n";
            }
            out << edges.size() << "\n";
            for (int i = 0; i < (int) edges.size(); ++i) {
                out << "e " << i << " " << edges[i].from << " " << edges[i].to << " " <<
                    edges[i].lib << " " << edges[i].weight << " " << edges[i].len;

                if (libs[edges[i].lib].type == Lib::RNA_PAIR || libs[edges[i].lib].type == Lib::RNA_SPLIT_30 ||
                        libs[edges[i].lib].type == Lib::RNA_SPLIT_50) {
                    out << " " << " \"coord: " <<
                        edges[i].coordBegin1 << " " << edges[i].coordEnd1 << " " <<
                        edges[i].coordBegin2 << " " << edges[i].coordEnd2;
                    if (edges[i].chr_name != "") {
                        out << " chr_name: " << edges[i].chr_name << " \"\n";
                    } else {
                        out << "\"\n";
                    }
                } else {
                    if (edges[i].info != "") {
                        out << " \"" << edges[i].info << "\"\n";
                    } else {
                        out << "\n";
                    }
                }
            }

            out.close();
        }

        ContigGraph ContigGraph::read(std::string fileName) {
            INFO("start read graph \"" << fileName << "\"");

            ContigGraph g;
            std::ifstream in(fileName);

            if (!in) {
                WARN("file " << fileName << " unable to open");
                return g;
            }

            size_t ln;
            size_t vn;
            size_t en;

            in >> ln;
            DEBUG("lib num=" << ln);

            g.libs.resize(ln);
            for (int i = 0; i < ln; ++i) {
                char c;
                int id;
                std::string type;
                in >> c >> id >> g.libs[i].color >> g.libs[i].name >> type;
                g.libs[i] = Lib(g.libs[i].color, g.libs[i].name, type);
            }

            in >> vn;
            DEBUG("vertex num=" << vn);
            g.graph.resize(vn);
            g.graphR.resize(vn);

            int mxT = 0;

            for (int i = 0; i < vn; ++i) {
                char c;
                in >> c;
                unsigned int v;
                in >> v;
                mxT = std::max((unsigned) mxT, v + 1);
                g.targets.resize(mxT);
                in >> g.targets[v].name;
                in >> g.targets[v].len;
                //TODO: read extra info to /dev/null ?
                g.targetId[g.targets[v].name] = v;
            }

            in >> en;
            DEBUG("edge num=" << en);
            g.edges.resize(en);
            std::string tmp;
            getline(in, tmp);

            for (int i = 0; i < en; ++i) {
                char c;
                int id;
                std::string curLine;
                getline(in, curLine);
                std::stringstream ss(curLine);
                ss >> c >> id >> g.edges[i].from >> g.edges[i].to >> g.edges[i].lib >> g.edges[i].weight >> g.edges[i].len;
                g.graph[g.edges[i].from].push_back(i);
                g.graphR[g.edges[i].to].push_back(i);
                //TODO: read extra info to /dev/null ?
                std::string tmp;
                if (ss >> tmp) {
                    ss >> g.edges[i].coordBegin1 >> g.edges[i].coordEnd1 >> g.edges[i].coordBegin2 >> g.edges[i].coordEnd2;
                }
                if (ss >> tmp) {
                    ss >> g.edges[i].chr_name;
                }
            }

            in.close();
            return g;
        }

        std::string ContigGraph::getLibColor(int l) {
            TRACE("getLibColor l=" << l << " : " << libs[l].color);
            return libs[l].color;
        }

        ContigGraph::Lib::Type ContigGraph::getLibType(int l) {
            TRACE("getLibType l=" << l << " : " << Lib::typeToStr[libs[l].type]);
            return libs[l].type;
        }

        void ContigGraph::setEdgeChr(int e, std::string name) {
            TRACE("setEdgeChr e=" << e << " name=" << name);
            edges[e].chr_name = name;
        }

        int ContigGraph::getTargetId(std::string name) {
            return targetId[name];
        }

        std::vector<ContigGraph::Edge> ContigGraph::getEdgesBetween(int v, int u) {
            std::vector<Edge> res;
            for (int id : graph[v]) {
                if (edges[id].to == u) {
                    res.push_back(edges[id]);
                }
            }
            if (res.size() > 500) {
                WARN("edeges betwwen = " << res.size());
                for (int i = 0; i < (int)res.size(); ++i) {
                    WARN("edge id" << res[i].id << " v=" << res[i].from << "(" <<res[i].coordBegin1 << ", "
                                   << res[i].coordEnd1 << ")"
                                   << "  u=" << res[i].to << "(" << res[i].coordBegin2 << ", " << res[i].coordEnd2 << ")");
                }
            }
            return res;
        }

        int ContigGraph::addEdge(int v1, int v2, std::pair<int, int> c1, std::pair<int, int> c2) {
            int e = (int) edges.size();
            edges.push_back(Edge(e, v1, v2, (int) libs.size() - 1, 0, c1.first, c1.second,
                                 c2.first, c2.second));

            graph[v1].push_back(e);
            graphR[v2].push_back(e);

            edges[e].weight = 1;
            edges[e].coordBegin1 = c1.first;
            edges[e].coordEnd1 = c1.second;
            edges[e].coordBegin2 = c2.first;
            edges[e].coordEnd2 = c2.second;

            return e;
        }

        void ContigGraph::incEdge(int e, std::pair<int, int> c1, std::pair<int, int> c2) {
            edges[e].weight += 1;
            edges[e].coordBegin1 = c1.first;
            edges[e].coordEnd1 = c1.second;
            edges[e].coordBegin2 = c2.first;
            edges[e].coordEnd2 = c2.second;
        }

        int ContigGraph::addEdge(int v1, int v2, double w, int len, std::string info) {
            int e = (int) edges.size();
            edges.push_back(Edge(e, v1, v2, (int) libs.size() - 1, w, 0, 0,
                                 0, 0));
            edges[edges.size() - 1].len = len;
            edges[edges.size() - 1].info = info;

            graph[v1].push_back(e);
            graphR[v2].push_back(e);
            return e;
        }
    }
}