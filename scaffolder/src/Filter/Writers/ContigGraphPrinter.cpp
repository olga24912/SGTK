//
// Created by olga on 23.10.16.
//

/*#include "ContigGraphPrinter.h"

void ContigGraphPrinter::writeFullGraphDotFormat(ContigGraph *g, string fileName) {
    cerr << "start write graph dot" << endl;
    ofstream out(fileName);

    out << "digraph {\n";

    for (int i = 0; i < (g->targetName).size(); ++i) {
        writeOneVertex(g, out, i);
    }

    for (int v = 0; v < (int)(g->graph).size(); ++v) {
        int vId = (g->idByV)[v];
        if (!isGoodVertex(g, v)) continue;
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

void ContigGraphPrinter::writeBigComponent(ContigGraph *g, int minSize, string fileName) {
    vector<int> drawV = vertexInBigComponents(g, minSize);
    if (drawV.size() < 2) return;

    writeThisVertex(g, drawV, fileName);
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

void ContigGraphPrinter::writeAlongPath(ContigGraph *g, int libId, int dist, int minSize, string fileName) {
    int n = (g->getVertexCount());
    int *col = new int[n];
    int cur = findComponent(g, col, IsGoodEdge(libId));

    vector< vector<int> > parts(cur);
    for (int i = 0; i < n; ++i) {
        if (col[i] == 0) continue;
        parts[col[i]].push_back(i);
    }

    for (int i = 1; i < cur; ++i) {
        if (parts[i].size() < minSize) continue;
        string fn = fileName;
        stringstream ss;
        ss << i;
        fn += ss.str();
        cerr << fn << endl;

        vector<int> res;
        for (int j = 0; j < (int)parts[i].size(); ++j) {
            vector<int> local = findAllVert(g, dist, parts[i][j]);
            for (int g = 0; g < (int)local.size(); ++g) {
                res.push_back(local[g]);
            }
        }

        sort(res.begin(), res.end());
        res.resize(unique(res.begin(), res.end()) - res.begin());

        writeThisVertex(g, res, fn);
    }
    delete col;
}
*/