#include <bits/stdc++.h>
#include "ContigGraph/ContigGraph.h"
#include "Logger/log_writers.hpp"

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

    create_console_logger("log.properties");

    INFO("start add info file \""<< infoFileName << "\" to graph \"" << grFileName
                                << "\" with new lib name \"" << libName << "\" and with new lib color " << color);

    ContigGraph graph = ContigGraph::read(grFileName);
    ifstream infoin(infoFileName);
    graph.newLib(libName, color, ContigGraph::Lib::Type::SCAFF);

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
                graph.incEdgeWeight(v, u, 0, 0, 0, 0);
                graph.incEdgeWeight(u ^ 1, v ^ 1, 0, 0, 0, 0);
            }

            v = u;
        }
    }

    graph.write("out.gr");

    infoin.close();
}
