//
// Created by olga on 09.10.16.
//

#include "Graph.h"

string Graph::genRandomColor() {
    int color[3] = {rand() % 256, rand() % 256, rand() % 256};
    string res = "#";
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
    res = "color = \"" + res + "\"";
    return res;
}

int Graph::getLibNum() {
    return libNum;
}
