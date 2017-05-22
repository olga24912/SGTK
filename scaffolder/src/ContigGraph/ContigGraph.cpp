#include <sstream>
#include "ContigGraph.h"

int ContigGraph::getLibNum() {
    return (int)libs.size();
}

int ContigGraph::getVertexCount() {
    return (int)graph.size();
}

int ContigGraph::getTargetLength(int v) const {
    return targets[v].len;
}

int ContigGraph::addVertex(int id, std::string name, int len) {
    assert(id == (int)graph.size());
    graph.push_back(std::vector<int>());
    graphR.push_back(std::vector<int>());
    targetId[name] = id;

    targets.resize((size_t)id + 1);
    targets[id].name = name;
    targets[id].len = len;

    vrtsToEdge.push_back(std::unordered_map<int, std::vector<int> >());
    return id;
}

int ContigGraph::incEdgeWeight(int v, int u) {
    assert(libs.size() > 0);
    assert(v < vrtsToEdge.size());
    int e;
    if (vrtsToEdge[v].count(u) == 0) {
        e = (int)edgeWeight.size();
        edgeWeight.push_back(0);
        edgeExtraInfo.push_back("");
        to.push_back(u);
        from.push_back(v);
        edgeLib.push_back((int)libName.size() - 1);

        vrtsToEdge[v][u] = e;

        graph[v].push_back(e);
        graphR[u].push_back(e);
    } else {
        e = vrtsToEdge[v][u];
    }

    edgeWeight[e] += 1;
    return e;
}

void ContigGraph::newLib(std::string name, std::string color) {
    libs.push_back(Lib(color, name));

    for (int i = 0; i < (int)vrtsToEdge.size(); ++i) {
        vrtsToEdge[i].clear();
    }
}

std::vector<int> ContigGraph::getEdges(int v) {
    return graph[v];
}

std::vector<int> ContigGraph::getEdgesR(int v) {
    return graphR[v];
}

int ContigGraph::getToVertex(int e) {
    return edges[e].to;
}

int ContigGraph::getFromVertex(int e) {
    return edges[e].from;
}

void ContigGraph::write(std::string fileName) {
    std::ofstream out(fileName);

    std::cerr << "Write -> gr" << " " << fileName << std::endl;
    std::cerr << edges.size() << std::endl;

    out << libs.size() << "\n";
    for (int i = 0; i < (int)libs.size(); ++i) {
        out << "l " << i << " " << libs[i].color << " " << libs[i].name << "\n";
    }
    out << graph.size() << "\n";
    for (int i = 0; i < (int)graph.size(); ++i) {
        out << "v " << i << " " << targets[i].name << " " << targets[i].len << "\n";
    }
    out << edges.size() << "\n";
    for (int i = 0; i < (int)edges.size(); ++i) {
        out << "e " << i << " " << edges[i].from << " " << edges[i].to << " " <<
            edges[i].lib << " " << edges[i].weight << " coord: " << edges[i].coordBegin1 << " " <<
            edges[i].coordEnd1 << " " << edges[i].coordBegin2 << " " << edges[i].coordEnd2 << "\n";
    }

    out.close();
}

ContigGraph ContigGraph::read(std::string fileName) {
    ContigGraph g;
    std::ifstream in(fileName);
    size_t ln;
    size_t vn;
    size_t en;
    in >> ln;
    std::cerr << ln << std::endl;
    g.libs.resize(ln);
    for (int i = 0; i < ln; ++i) {
        char c;
        int id;
        in >> c >> id >> g.libs[i].color >> g.libs[i].name;
    }

    in >> vn;
    std::cerr << vn << std::endl;
    g.graph.resize(vn);
    g.graphR.resize(vn);
    g.vrtsToEdge.resize(vn);

    int mxT = 0;

    for (int i = 0; i < vn; ++i) {
        char c;
        in >> c;
        unsigned int v;
        in >> v;
        mxT = std::max((unsigned)mxT, v + 1);
        g.targets.resize(mxT);
        in >> g.targets[v].name;
        in >> g.targets[v].len;
        g.targetId[g.targets[v].name] = v;
    }

    in >> en;
    std::cerr << en << std::endl;
    g.edges.resize(en);

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
        ss >> tmp >> g.edges[i].coordBegin1 >> g.edges[i].coordEnd1 >> g.edges[i].coordBegin2 >> g.edges[i].coordEnd2;
    }

    in.close();
    return g;
}

std::string ContigGraph::getTargetName(int v) {
    return targets[v].name;
}

int ContigGraph::getEdgeWeight(int e) {
    return edges[e].weight;
}

int ContigGraph::getEdgeLib(int e) {
    return edges[e].lib;
}

std::string ContigGraph::getLibColor(int l) {
    return libs[l].color;
}

std::string ContigGraph::getLibName(int l) {
    return libs[l].name;
}

int ContigGraph::getTargetId(std::string name) {
    return targetId[name];
}

std::string ContigGraph::getEdgeInfo(int e) {
    std::stringstream ss;
    ss << "coord: " << edges[e].coordBegin1 << "-" << edges[e].coordEnd1 << "\n" << edges[e].coordBegin2 << "-" << edges[e].coordEnd2;
    return ss.str();
}
