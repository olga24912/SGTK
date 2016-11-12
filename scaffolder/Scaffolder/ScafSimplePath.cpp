#include "ScafSimplePath.h"
#include "../Scaffolder/ScafSimplePath.h"

void ScafSimplePath::evaluate(ContigGraph* graph, string contigFile, string out) {
    this->graph = graph;
    this->contigFile = contigFile;

    cerr << "save" << endl;
    saveContigs();
    cerr << "find path" << endl;
    findPaths();
    cerr << "write new contig" << endl;
    writeNewContigs(out);
}

void ScafSimplePath::saveContigs() {
    SeqFileIn seqFileIn(contigFile.c_str());
    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;

    readRecords(ids, seqs, seqFileIn);

    for (unsigned i = 0; i < length(ids); ++i) {
        string seq = SeqanUtils::dna5ToString(toCString(seqs[i]), length(seqs[i]));
        contigs.push_back(seq);
        contigs.push_back(createRevCompl(seq));
    }
}

void ScafSimplePath::findPaths() {
    int n = graph->getVertexCount();
    vector<int> used(n, 0);
    next.resize(n, -1);

    for (int i = 0; i < n; ++i) {
        if (used[i] == 0 && graph->isGoodVertex(i)) {
            dfsPath(i, next, used);
        }
        cerr << used[i] << " " << next[i] << endl;
    }
}

void ScafSimplePath::writeNewContigs(string fileName) {
    SeqFileOut out(fileName.c_str());

    StringSet<seqan::CharString> ids;
    StringSet<seqan::CharString> seqs;

    int n = graph->getVertexCount();
    vector<int> used(n, 0);
    vector<int> first(n, 1);
    for (int i = 0; i < n; ++i) {
        if (next[i] != -1) {
            first[next[i]] = 0;
        }
    }

    for (int i = 0; i < n; ++i) {
        cerr << i << endl;
        if (first[i] && used[i] == 0) {
            string seq = contigs[i];
            used[i] = 1;
            used[i ^ 1] = 1;

            int x = next[i];
            while (x != -1) {
                for (int j = 0; j < GAP_SIZE; ++j) {
                    seq += 'N';
                }
                seq += contigs[x];
                used[x] = 1;
                used[x^1] = 1;
                x = next[x];
            }

            stringstream ss;
            ss << i;
            string readName = "path" + ss.str();
            cerr << readName << endl;
            appendValue(ids, CharString(readName.c_str()));
            cerr << seq << " " << seq.size() <<  endl;
            appendValue(seqs, Dna5String(seq));
        }
    }

    writeRecords(out, ids, seqs);
}

string ScafSimplePath::createRevCompl(string s) {
    string res = s;
    reverse(res.begin(), res.end());
    for (int i = 0; i < (int)res.size(); ++i) {
        if (res[i] == 'A') {
            res[i] = 'T';
        } else if (res[i] == 'T') {
            res[i] = 'A';
        } else if (res[i] == 'C') {
            res[i] = 'G';
        } else if (res[i] == 'G') {
            res[i] = 'C';
        }
    }
    return res;
}

void ScafSimplePath::dfsPath(int v, vector<int> &next, vector<int> &used) {
    used[v] = 1;

    if (!graph->isGoodVertex(v)) {
        used[v] = 3;
        return;
    }

    int tu = -1;
    vector<int> edges = graph->getEdges(v);
    for (int i = 0; i < (int)edges.size(); ++i) {
        int e = edges[i];
        int u = graph->getToVertex(e);

        if (!graph->isGoodEdge(e) || !graph->isGoodVertex(u)) continue;
        if (tu == -1) {
            tu = u;
        }
        if (tu != u || used[u] == 1) {
            used[v] = 3;
            return;
        }
    }

    vector<int> edgeR = graph->getEdgesR(v);
    int fu = -1;

    for (int i = 0; i < (int)edgeR.size(); ++i) {
        int e = edgeR[i];
        int u = graph->getFromVertex(e);

        if (!graph->isGoodEdge(e) || !graph->isGoodVertex(u)) continue;
        if (fu == -1) {
            fu = u;
        }
        if (fu != u) {
            used[v] = 3;
            return;
        }
    }

    next[v] = tu;
    if (tu != -1 && used[tu] == 0) dfsPath(tu, next, used);

    if (tu != -1 && used[tu] == 3) {
        next[v] = -1;
    }

    used[v] = 2;
}
