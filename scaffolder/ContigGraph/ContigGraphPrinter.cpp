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

    writeThisVertex(g, drawV, fileName);
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
    int cntEdge = 0;
    for (int i = 0; i < (g->graph)[v].size(); ++i) {
        int e = (g->graph)[v][i];
        int uid = g->idByV[(g->to)[e]];
        if ((g->edgeWeight)[e] >= (g->libMinEdgeWight)[g->edgeLib[e]] &&
                (g->targetLen)[uid] >= g->minContigLen) {
            ++cntEdge;
        }
    }

    for (int i = 0; i < (g->graphR)[v].size(); ++i) {
        int e = (g->graphR)[v][i];
        int uid = g->idByV[(g->from)[e]];
        if ((g->edgeWeight)[e] >= (g->libMinEdgeWight)[g->edgeLib[e]] &&
                (g->targetLen)[uid] >= g->minContigLen) {
            ++cntEdge;
        }
    }

    if (cntEdge == 0) return;

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

void ContigGraphPrinter::writeBigComponent(ContigGraph *g, int minSize, string fileName) {
    vector<int> drawV = vertexInBigComponents(g, minSize);
    if (drawV.size() < 2) return;

    writeThisVertex(g, drawV, fileName);
}

void ContigGraphPrinter::writeThisVertex(ContigGraph *g, vector<int> &drawV, string fileName) {
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

vector<int> ContigGraphPrinter::vertexInBigComponents(ContigGraph *g, int size) {
    vector<int> res;

    int n = (g->getVertexCount());
    int col[n];
    for (int i = 0; i < n; ++i) {
        col[i] = 0;
    }
    int cur = 1;

    for (int i = 0; i < n; ++i) {
        if (col[i] == 0) {
            dfsFindComponent(g, col, cur, i);
            ++cur;
        }
    }

    vector<int> cntCol(cur, 0);

    for (int i = 0; i < n; ++i) {
        cntCol[col[i]]++;
    }

    for (int  i  = 0; i < n; ++i) {
        if (cntCol[i] >= size) {
            res.push_back(i);
        }
    }

    return res;
}

void ContigGraphPrinter::dfsFindComponent(ContigGraph *g, int *color, int currentCol, int v) {
    color[v] = currentCol;

    for (int i = 0; i < (int)g->graph[v].size(); ++i) {
        int e = g->graph[v][i];
        int u = g->to[e];
        if (color[u] == 0) {
            dfsFindComponent(g, color, currentCol, u);
        }
    }

    for (int i = 0; i < (int)g->graphR[v].size(); ++i) {
        int e = g->graphR[v][i];
        int u = g->from[e];
        if (color[u] == 0) {
            dfsFindComponent(g, color, currentCol, u);
        }
    }
}
