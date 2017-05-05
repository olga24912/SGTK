#include <sstream>
#include "ContigGraph.h"

int ContigGraph::getLibNum() {
    return (int)libName.size();
}

int ContigGraph::getVertexCount() {
    return (int)graph.size();
}

int ContigGraph::getTargetLength(int v) const {
    return targetLen[v];
}

int ContigGraph::addVertex(int id, std::string name, int len) {
    assert(id == (int)graph.size());

    graph.push_back(std::vector<int>());
    graphR.push_back(std::vector<int>());

    targetId[name] = id;

    targetName.resize((size_t)id + 1);
    targetName[id] = name;

    targetLen.resize((size_t)id + 1);
    targetLen[id] = len;

    vrtsToEdge.push_back(std::unordered_map<int, int>());
    return id;
}

int ContigGraph::incEdgeWeight(int v, int u) {
    assert(libName.size() > 0);
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
    libColor.push_back(color);
    libName.push_back(name);

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
    return to[e];
}

int ContigGraph::getFromVertex(int e) {
    return from[e];
}

void ContigGraph::write(std::string fileName) {
    std::ofstream out(fileName);

    std::cerr << "Write -> gr" << " " << fileName << std::endl;
    std::cerr << edgeWeight.size() << std::endl;

    out << libName.size() << "\n";
    for (int i = 0; i < (int)libName.size(); ++i) {
        out << "l " << i << " " << libColor[i] << " " << libName[i] << "\n";
    }
    out << graph.size() << "\n";
    for (int i = 0; i < (int)graph.size(); ++i) {
        out << "v " << i << " " << targetName[i] << " " << targetLen[i] << "\n";
    }
    out << edgeWeight.size() << "\n";
    for (int i = 0; i < (int)edgeWeight.size(); ++i) {
        out << "e " << i << " " << from[i] << " " << to[i] << " " <<
            edgeLib[i] << " " << edgeWeight[i] << " " << edgeExtraInfo[i] << "\n";
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
    g.libName.resize(ln);
    g.libColor.resize(ln);
    for (int i = 0; i < ln; ++i) {
        char c;
        int id;
        in >> c >> id >> g.libColor[i] >> g.libName[i];
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
        g.targetName.resize(mxT);
        in >> g.targetName[v];
        g.targetLen.resize(mxT);
        in >> g.targetLen[v];
        g.targetId[g.targetName[v]] = v;
    }

    in >> en;
    std::cerr << en << std::endl;
    g.to.resize(en);
    g.from.resize(en);
    g.edgeLib.resize(en);
    g.edgeWeight.resize(en);
    g.edgeExtraInfo.resize(en);

    for (int i = 0; i < en; ++i) {
        char c;
        int id;
        std::string curLine;
        getline(in, curLine);
        std::stringstream ss(curLine);
        ss >> c >> id >> g.from[i] >> g.to[i] >> g.edgeLib[i] >> g.edgeWeight[i];
        g.graph[g.from[i]].push_back(i);
        g.graphR[g.to[i]].push_back(i);
        getline(ss, g.edgeExtraInfo[id]);
    }

    in.close();
    return g;
}

std::string ContigGraph::getTargetName(int v) {
    return targetName[v];
}

int ContigGraph::getEdgeWeight(int e) {
    return edgeWeight[e];
}

int ContigGraph::getEdgeLib(int e) {
    return edgeLib[e];
}

std::string ContigGraph::getLibColor(int l) {
    return libColor[l];
}

std::string ContigGraph::getLibName(int l) {
    return libName[l];
}

int ContigGraph::getTargetId(std::string name) {
    return targetId[name];
}

void ContigGraph::setEdgeInfo(int e, std::string info) {
    edgeExtraInfo[e] = info;
}

std::string ContigGraph::getEdgeInfo(int e) {
    return edgeExtraInfo[e];
}
