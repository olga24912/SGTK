#include <sstream>
#include "ContigGraph.h"
#include "cmath"

namespace contig_graph {
    const std::string ContigGraph::Lib::typeToStr[] = {"REF", "DNA_PAIR", "RNA_PAIR", "RNA_SPLIT_50", "RNA_SPLIT_30",
                                                       "SCAFF"};

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

        vrtsToEdge.push_back(std::unordered_map<int, std::vector<int> >());
        return id;
    }

    int ContigGraph::incEdgeWeight(int v, int u, int cb1, int ce1, int cb2, int ce2) {
        assert(libs.size() > 0);
        assert(v < vrtsToEdge.size());

        DEBUG("increment edge Weight between v=" << v << "(" << cb1 << "-" << ce1 << ") u=" << u << "(" << cb2 << "-"
                                                 << ce2 << ")");

        int e = -1;
        std::vector<int> candidates = vrtsToEdge[v][u];

        for (int ec : candidates) {
            if (std::fabs(cb1 - edges[ec].coordBegin1) < maxClusterSize &&
                std::fabs(cb2 - edges[ec].coordBegin2) < maxClusterSize) {
                e = ec;
            }
        }

        if (e == -1) {
            e = (int) edges.size();
            edges.push_back(Edge(e, v, u, (int) libs.size() - 1, 0, cb1, ce1, cb2, ce2));

            vrtsToEdge[v][u].push_back(e);

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

        for (int i = 0; i < (int) vrtsToEdge.size(); ++i) {
            vrtsToEdge[i].clear();
        }
    }

    std::vector<int> ContigGraph::getEdges(int v) {
        TRACE("getEdges v=" << v);
        return graph[v];
    }

    std::vector<int> ContigGraph::getEdgesR(int v) {
        TRACE("getEdgesR v=" << v);
        return graphR[v];
    }

    int ContigGraph::getToVertex(int e) {
        TRACE("get v (e(" << e << "):u(" << edges[e].from << ")->V(" << edges[e].to << "))");
        return edges[e].to;
    }

    int ContigGraph::getFromVertex(int e) {
        TRACE("get U (e(" << e << "):U(" << edges[e].from << ")->v(" << edges[e].to << "))");
        return edges[e].from;
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
                edges[i].lib << " " << edges[i].weight << " coord: " << edges[i].coordBegin1 << " " <<
                edges[i].coordEnd1 << " " << edges[i].coordBegin2 << " " << edges[i].coordEnd2;
            if (edges[i].chr_name != "") {
                out << " chr_name: " << edges[i].chr_name << "\n";
            } else {
                out << "\n";
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
        g.vrtsToEdge.resize(vn);

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
            ss >> c >> id >> g.edges[i].from >> g.edges[i].to >> g.edges[i].lib >> g.edges[i].weight;
            g.graph[g.edges[i].from].push_back(i);
            g.graphR[g.edges[i].to].push_back(i);
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

    std::string ContigGraph::getTargetName(int v) {
        TRACE("getTargetName v=" << v << " : " << targets[v].name);
        return targets[v].name;
    }

    int ContigGraph::getEdgeWeight(int e) {
        TRACE("getEdgeWeight e=" << e << " : " << edges[e].weight);
        return edges[e].weight;
    }

    int ContigGraph::getEdgeLib(int e) {
        TRACE("getEdgeLib e=" << e << " : " << edges[e].lib);
        return edges[e].lib;
    }

    std::string ContigGraph::getLibColor(int l) {
        TRACE("getLibColor l=" << l << " : " << libs[l].color);
        return libs[l].color;
    }

    std::string ContigGraph::getLibName(int l) {
        TRACE("getLibName l=" << l << " : " << libs[l].name);
        return libs[l].name;
    }

    int ContigGraph::getTargetId(std::string name) {
        TRACE("getTargetId name=" << name << " : " << targetId[name]);
        return targetId[name];
    }

    std::string ContigGraph::getEdgeInfo(int e) {
        std::stringstream ss;
        ss << edges[e].chr_name << "\n";
        ss << "coord: " << edges[e].coordBegin1 << "-" << edges[e].coordEnd1 << "\n" << edges[e].coordBegin2 << "-"
           << edges[e].coordEnd2;
        TRACE("getEdgeInfo e=" << e << " : " << ss.str());
        return ss.str();
    }

    ContigGraph::Lib::Type ContigGraph::getLibType(int l) {
        TRACE("getLibType l=" << l << " : " << Lib::typeToStr[libs[l].type]);
        return libs[l].type;
    }

    int ContigGraph::getEdgeCoordB1(int e) {
        TRACE("getEdgeCoordB1 e=" << e << " : " << edges[e].coordBegin1);
        return edges[e].coordBegin1;
    }

    int ContigGraph::getEdgeCoordE1(int e) {
        TRACE("getEdgeCoordE1 e=" << e << " : " << edges[e].coordEnd1);
        return edges[e].coordEnd1;
    }

    int ContigGraph::getEdgeCoordB2(int e) {
        TRACE("getEdgeCoordB2 e=" << e << " : " << edges[e].coordBegin2);
        return edges[e].coordBegin2;
    }

    int ContigGraph::getEdgeCoordE2(int e) {
        TRACE("getEdgeCoordE2 e=" << e << " : " << edges[e].coordEnd2);
        return edges[e].coordEnd2;
    }

    void ContigGraph::setEdgeChr(int e, std::string name) {
        TRACE("setEdgeChr e=" << e << " name=" << name);
        edges[e].chr_name = name;
    }
}
