//
// Created by olga on 23.10.16.
//

#include "ContigGraphPrinter.h"

void ContigGraphPrinter::writeFullGraphDotFormat(ContigGraph *g, string fileName) {
    cerr << "start write graph dot" << endl;
    ofstream out(fileName);

    out << "digraph {\n";

    for (int i = 0; i < (g->targetName).size(); ++i) {
        writeOneVertex(g, out, i);
    }

    for (int v = 0; v < (int)(g->graph).size(); ++v) {
        int vId = (g->idByV)[v];
        if ((g->targetLen)[vId] < (g->minContigLen)) continue;
        for (int j = 0; j < (int)(g->graph)[v].size(); ++j) {
            int e = (g->graph)[v][j];
            writeOneEdge(g, out, v, e);
        }
    }
    out << "}\n";

    out.close();
}

void ContigGraphPrinter::writeAllLocalGraphDotFormat(ContigGraph *g, int dist) {
    for (int i = 0; i < (g->targetName).size(); ++i) {
        string name = "";
        int x = i;
        while (x > 0) {
            name += '0' + (x%10);
            x /= 10;
        }
        reverse(name.begin(), name.end());
        name = "g" + name;
        writeLocalGraph(g, dist, i, name);
    }
}


void ContigGraphPrinter::writeLocalGraph(ContigGraph *g, int dist, int v, string fileName) {
    vector<int> drawV = findAllVert(g, dist, v);
    if (drawV.size() < 2) return;

    ofstream out(fileName);
    out << "digraph {\n";
    for (int i = 0; i < (int)drawV.size(); ++i) {
        writeOneVertex(g, out, drawV[i]);
    }

    for (int i = 0; i < (int)drawV.size(); ++i) {
        int v = drawV[i];
        for (int j = 0; j < (g->graph)[v].size(); ++j) {
            int e = (g->graph)[v][j];
            int u = (g->to)[e];

            int was = 0;
            for (int h = 0; h < (int)drawV.size(); ++h) {
                if (drawV[h] == u) was = 1;
            }
            if (was) {
                writeOneEdge(g, out, v, e);
            }
        }
    }

    out << "}\n";
    out.close();
}

void ContigGraphPrinter::writeOneEdge(ContigGraph *g, ofstream &out, int v, int e) {
    int vId = (g->idByV)[v];
    int u = (g->to)[e];
    int uId = (g->idByV)[u];
    if ((g->targetLen)[uId] < (g->minContigLen)) return;
    if ((g->edgeWeight)[e] < (g->libMinEdgeWight)[(g->edgeLib[e])]) return;
    out << "    \"" << (g->targetName)[vId] << "\" -> \"" << (g->targetName)[uId] << "\" [ ";
    out << "color = \"" << (g->libColor)[(g->edgeLib)[e]] << "\", ";
    out << "penwidth = "<< 1 + (int)log10((g->edgeWeight)[e]) << ", ";
    out << "label = " << "\"" << (g->libName)[(g->edgeLib)[e]] << "\n weight = " << (g->edgeWeight)[e] << "\" ]\n";
}


void ContigGraphPrinter::writeOneVertex(ContigGraph *g, ofstream &out, int v) {
    int vId = (g->idByV)[v];
    if ((g->targetLen)[vId] < (g->minContigLen)) return;
    out << "    \"" << (g->targetName)[vId] << "\"[label=\" " << (g->targetName)[vId] <<
    " id = " << v << "\nlen = "
    << (g->targetLen)[vId]
    << ", cover = "<< (g->targetCoverage)[vId] <<"\"];\n";
}

vector<int> ContigGraphPrinter::findAllVert(ContigGraph *g, int dist, int v) {
    vector<int> res;
    res.push_back(v);
    if ((g->targetLen)[(g->idByV)[v]] < (g->minContigLen)) return res;
    if (dist == 0) return res;
    for (int i = 0; i < (g->graph)[v].size(); ++i) {
        int e = (g->graph)[v][i];
        if ((g->edgeWeight)[e] < (g->libMinEdgeWight)[(g->edgeLib)[e]]) continue;
        int u = (g->to)[e];
        vector<int> add = findAllVert(g, dist - 1, u);
        for (int j = 0; j < (int)add.size(); ++j) {
            res.push_back(add[j]);
        }
    }

    for (int i = 0; i < (g->graphR)[v].size(); ++i) {
        int e = (g->graphR)[v][i];
        if ((g->edgeWeight)[e] < (g->libMinEdgeWight)[(g->edgeLib)[e]]) continue;
        int u = (g->from)[e];
        vector<int> add = findAllVert(g, dist - 1, u);
        for (int j = 0; j < (int)add.size(); ++j) {
            res.push_back(add[j]);
        }
    }

    sort(res.begin(), res.end());
    res.resize(unique(res.begin(), res.end()) - res.begin());

    return res;
}
