#include "GraphBuilder.h"

void GraphBuilder::setGraph(ContigGraph *graph) {
    DEBUG("setGraph");
    GraphBuilder::graph = graph;
    graph->newLib(libName, getLibColor(), getLibType());
}

void GraphBuilder::setLibName(std::string libName, std::string path) {
    DEBUG("setLibName libName=" << libName << " path=" << path);
    this->libName = libName;
    GraphBuilder::path = path + "/" + libName;

    std::string command = "mkdir " + GraphBuilder::path;
    system(command.c_str());
    DEBUG("command: " << command);
}

std::string GraphBuilder::colorToString(int *color) {
    std::string res = "#";
    for (int i = 0; i < 3; ++i) {
        if (color[i] / 16 < 10) {
            res += (color[i] / 16) + '0';
        } else {
            res += (color[i] / 16) - 10 + 'a';
        }

        if (color[i] % 16 < 10) {
            res += (color[i] % 16) + '0';
        } else {
            res += (color[i] % 16) - 10 + 'a';
        }
    }
    TRACE("colorToString color=(" << color[0] << " " << color[1] << " " << color[2] << ") : " << res);
    return res;
}

void GraphBuilder::setSamFileWriter() {
    DEBUG("setSamFileWriter");
    this->samFileWriter = SamFileWriteEdge(path);
}
