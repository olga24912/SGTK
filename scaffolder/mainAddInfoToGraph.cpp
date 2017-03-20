#include <bits/stdc++.h>
#include "ContigGraph/ContigGraph.h"

using namespace std;

//add edges form info to gr
//argv[1] - name of info file
//argv[2] - name of gr file
int main(int argv, char** argc) {
    string infoFileName = argc[1];
    string grFileName = argc[2];

    ContigGraph graph = ContigGraph::read(grFileName);

    ifstream infoin(infoFileName);

    graph.newLib("infoLib", "#ffae00");

    string cur;
    while (getline(infoin, cur)) {
        stringstream ss(cur);

        string x;
        ss >> x;
        string contigName;
        int id;
        char dir;
        char tr;

        int v = -1;
        while (ss >> contigName >> id >> dir >> tr) {
            contigName = contigName.substr(1);
            int u = graph.getTargetId(contigName);
            if (dir == '-') {
                u ^= 1;
            }

            if (v != -1) {
                graph.incEdgeWeight(v, u);
                graph.incEdgeWeight(u ^ 1, v ^ 1);
            }

            v = u;
        }
    }

    graph.write("out.gr");

    infoin.close();
}
