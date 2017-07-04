#include <bits/stdc++.h>
#include <Logger/log_writers.hpp>
#include "ContigGraph/ContigGraph.h"

using namespace std;

void create_console_logger(const string& log_props_file) {
    using namespace logging;

    //string log_props_file = cfg::get().log_filename;

    //if (!path::FileExists(log_props_file))
    //    log_props_file = path::append_path(dir, cfg::get().log_filename);

    logger *lg = create_logger(path::FileExists(log_props_file) ? log_props_file : "");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}
//add edges form both.path to gr
//argc[1] - name of both.path file
//argc[2] - name of gr file
//argc[3] - lib name
//argc[4] - color
int main(int argv, char** argc) {
    create_console_logger("../log.properties");
    string fileName = argc[1];
    string grFileName = argc[2];
    string libName = argc[3];
    string color = argc[4];

    ContigGraph graph = ContigGraph::read(grFileName);

    ifstream infoin(fileName);

    graph.newLib(libName, color, ContigGraph::Lib::Type::SCAFF);

    string cur;
    while (getline(infoin, cur)) {
        std::vector<int> vid;
        string curName = "";
        int i = 0;
        while (i <= (int)cur.size()) {
            if (i < cur.size() && cur[i] != '-' && cur[i] != '>') {
                curName += cur[i];
            } else if (cur[i] == '-' || i == (int)cur.size()){
                if (curName[curName.size() - 1] != ')') {
                    if (curName[curName.size() - 2] == '/') {
                        curName.resize(curName.size() - 2);
                        int v = graph.getTargetId(curName);
                        vid.push_back(v^1);
                    } else {
                        vid.push_back(graph.getTargetId(curName));
                    }
                }
                curName = "";
            }
            ++i;
        }
        for (int j = 1; j < (int)vid.size(); ++j) {
                graph.incEdgeWeight(vid[j - 1], vid[j], 0, 0, 0, 0);
                graph.incEdgeWeight(vid[j] ^ 1, vid[j - 1] ^ 1, 0, 0, 0, 0);
        }
    }

    graph.write("out.gr");

    infoin.close();
}
