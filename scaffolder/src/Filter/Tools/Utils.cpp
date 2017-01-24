//
// Created by olga on 04.12.16.
//

#include "Utils.h"

bool Utils::isInt(string s) {
    for (int i = 0; i < (int)s.size(); ++i) {
        if (!(s[i] >= '0' && s[i] <= '9')) {
            return false;
        }
    }
    return true;
}

int Utils::strToInt(string s) {
    stringstream ss(s);
    int c;
    ss >> c;
    return c;
}
