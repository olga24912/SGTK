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
        name = "gl" + name;
        writeLocalGraph(g, dist, i, name);
    }
}


void ContigGraphPrinter::writeLocalGraph(ContigGraph *g, int dist, int v, string fileName) {
    vector<int> drawV = findAllVert(g, dist, v);
    if (drawV.size() < 2) return;

    //cerr << drawV.size() << endl;

    writeThisVertex(g, drawV, fileName);
}

void ContigGraphPrinter::writeLocalSegGraph(ContigGraph *g, int dist, int vb, int ve, string fileName) {
    vector<int> drawV;
    for (int i = vb; i <= ve; ++i) {
        vector<int> lv = findAllVert(g, dist, i);
        for (int j = 0; j < lv.size(); ++j) {
            drawV.push_back(lv[j]);
        }
    }

    sort(drawV.begin(), drawV.end());

//    cerr << drawV.size() << endl;

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
    out << "label = " << "\"" << (g->libName)[(g->edgeLib)[e]] << "\n weight = " << (g->edgeWeight)[e];
    out << "\n id = "<< e << "\" ]\n";
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
    if ((g->targetLen)[(g->idByV)[v]] < (g->minContigLen)) return res;
    if (dist == 0) return res;
    res.push_back(v);
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
    vector<pair<int,int> > wieghtEdge;

    for (int i = 0; i < (int)drawV.size(); ++i) {
        int v = drawV[i];
        for (int j = 0; j < (g->graph)[v].size(); ++j) {
            int e = (g->graph)[v][j];
            int u = (g->to)[e];

            if ((g->edgeWeight)[e] < (g->libMinEdgeWight)[(g->edgeLib[e])]) continue;

            int was = 0;
            for (int h = 0; h < (int)drawV.size(); ++h) {
                if (drawV[h] == u) was = 1;
            }
            if (was) {
                wieghtEdge.push_back(make_pair(g->edgeWeight[e], e));
            }
        }
    }

    if (drawV.size() == 5 && wieghtEdge.size() == 4) return;
    if (fileName == "gl10012") {
        for (int i  = 0; i < drawV.size(); ++i) {
            cerr << drawV[i] << " ";
        }
        cerr << endl;
        for (int i = 0; i < wieghtEdge.size(); ++i) {
            cerr << wieghtEdge[i].first << " " << wieghtEdge[i].second << endl;
        }
        cerr << endl;
    }
    ofstream out(fileName);
    out << "digraph {\n";
    for (int i = 0; i < (int)drawV.size(); ++i) {
        writeOneVertex(g, out, drawV[i]);
    }

    sort(wieghtEdge.rbegin(), wieghtEdge.rend());
    for(int i = 0; i < (int)wieghtEdge.size(); ++i) {
        writeOneEdge(g, out, (g->from)[wieghtEdge[i].second], wieghtEdge[i].second);
    }

    out << "}\n";
    out.close();
}

vector<int> ContigGraphPrinter::vertexInBigComponents(ContigGraph *g, int size) {
    vector<int> res;
    int n = (g->getVertexCount());
    int *col = new int[n];
    int cur = findComponent(g, col);

    vector<int> cntCol(cur, 0);
    for (int i = 0; i < n; ++i) {
        cntCol[col[i]]++;
    }

    for (int  i  = 0; i < n; ++i) {
        if (col[i] != 0 && cntCol[col[i]] >= size) {
            res.push_back(i);
        }
    }

    delete col;
    return res;
}

int ContigGraphPrinter::findComponent(ContigGraph *g, int *col) {
    int n = (g->getVertexCount());
    int cur = 1;
    for (int i = 0; i < n; ++i) {
        col[i] = 0;
    }
    for (int i = 0; i < n; ++i) {
        if (col[i] == 0) {
            dfsFindComponent(g, col, cur, i);
            ++cur;
        }
    }
    return cur;
}

void ContigGraphPrinter::dfsFindComponent(ContigGraph *g, int *color, int currentCol, int v) {
    if (g->targetLen[g->idByV[v]] < g->minContigLen) return;
    color[v] = currentCol;

    for (int i = 0; i < (int)g->graph[v].size(); ++i) {
        int e = g->graph[v][i];
        int u = g->to[e];
        if (g->edgeWeight[e] < g->libMinEdgeWight[g->edgeLib[e]]) continue;
        if (color[u] == 0) {
            dfsFindComponent(g, color, currentCol, u);
        }
    }

    for (int i = 0; i < (int)g->graphR[v].size(); ++i) {
        int e = g->graphR[v][i];
        int u = g->from[e];
        if (g->edgeWeight[e] < g->libMinEdgeWight[g->edgeLib[e]]) continue;
        if (color[u] == 0) {
            dfsFindComponent(g, color, currentCol, u);
        }
    }
}

void ContigGraphPrinter::writeSplitBigComponent(ContigGraph *g, int minSize, string fileName) {
    int n = (g->getVertexCount());
    int *col = new int[n];
    int cur = findComponent(g, col);

    vector< vector<int> > parts(cur);
    for (int i = 0; i < n; ++i) {
        if (col[i] == 0) continue;
        parts[col[i]].push_back(i);
    }

    for (int i = 1; i < cur; ++i) {
        if (parts[i].size() >= minSize) {
            string fn = fileName;
            stringstream ss;
            ss << i;
            fn += ss.str();
            cerr << fn << endl;
            writeThisVertex(g, parts[i], fn);
        }
    }
    delete col;
}

