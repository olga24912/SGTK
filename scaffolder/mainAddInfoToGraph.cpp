#include <bits/stdc++.h>
#include "ContigGraph/ContigGraph.h"

using namespace std;

//add edges form info to gr
//argc[1] - name of info file
//argc[2] - name of gr file
//argc[3] - lib name
//argc[4] - color
int main(int argv, char** argc) {
    string infoFileName = argc[1];
    string grFileName = argc[2];
    string libName = argc[3];
    string color = argc[4];

    ContigGraph graph = ContigGraph::read(grFileName);

    ifstream infoin(infoFileName);

    graph.newLib(libName, color);

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
