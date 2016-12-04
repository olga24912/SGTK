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
        out << "v " << i << " " << (g->idByV)[i] << " " << (g->targetName)[(g->idByV)[i]] << " "
        << (g->targetCoverage)[(g->idByV)[i]] << " " << (g->targetLen)[(g->idByV)[i]] << "\n";
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
    g.libMinEdgeWight.resize(ln, 0);
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
    g.start.resize(vn);
    g.idByV.resize(vn);
    int mxT = 0;

    for (int i = 0; i < vn; ++i) {
        char c;
        in >> c;
        int v, vId;
        in >> v >> vId;
        g.idByV[v] = vId;
        mxT = max(mxT, vId + 1);
        g.vById.resize(mxT);
        g.vById[vId] = v;
        g.targetName.resize(mxT);
        in >> g.targetName[vId];
        g.targetCoverage.resize(mxT);
        in >> g.targetCoverage[vId];
        g.targetLen.resize(mxT);
        in >> g.targetLen[vId];
        g.targetId[g.targetName[vId]] = vId;
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
