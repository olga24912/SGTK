//
// Created by olga on 23.10.16.
//

#include "Serialization.h"

void Serialization::write(ContigGraph *g, string fileName) {
    ofstream out(fileName);
    cerr << "Write -> gr" << endl;

    out << (g->libName).size() << "\n";
    for (int i = 0; i < (int)(g->libName).size(); ++i) {
        out << "l " << i << " " << (g->libColor)[i] << " " << (g->libName)[i] << "\n";
    }
    out << (g->graph).size() << "\n";
    for (int i = 0; i < (int)(g->graph).size(); ++i) {
        out << "v " << i << " " << (g->targetName)[i] << " " << (g->targetLen)[i] << "\n";
    }
    out << (g->edgeWeight).size() << "\n";
    for (int i = 0; i < (int)(g->edgeWeight).size(); ++i) {
        out << "e " << i << " " << (g->from)[i] << " " << (g->to)[i] <<
                " " << (g->edgeLib)[i] << " " << (g->edgeWeight)[i] << "\n";
    }

    out.close();
}

ContigGraph Serialization::read(string fileName) {
    ContigGraph g;

    ifstream in(fileName);

    int ln;
    in >> ln;
    cerr << ln << endl;
    g.libName.resize(ln);
    for (int i = 0; i < ln; ++i) {
        char c;
        int id;
        string color;
        in >> c >> id >> color >> g.libName[i];
        g.libColor.push_back(color);
    }

    int vn;
    in >> vn;
    cerr << vn << endl;
    g.graph.resize(vn);
    g.graphR.resize(vn);

    int mxT = 0;

    for (int i = 0; i < vn; ++i) {
        char c;
        in >> c;
        int v;
        in >> v;
        mxT = max(mxT, v + 1);
        g.targetName.resize(mxT);
        in >> g.targetName[v];
        g.targetLen.resize(mxT);
        in >> g.targetLen[v];
        g.targetId[g.targetName[v]] = v;
    }

    int en;
    in >> en;
    cerr << en << endl;
    g.to.resize(en);
    g.from.resize(en);
    g.edgeLib.resize(en);
    g.edgeWeight.resize(en);

    for (int i = 0; i < en; ++i) {
        char c;
        int id;
        in >> c >> id >> g.from[i] >> g.to[i] >> g.edgeLib[i] >> g.edgeWeight[i];
        g.graph[g.from[i]].push_back(i);
        g.graphR[g.to[i]].push_back(i);
    }

    in.close();
    return g;
}
