//
// Created by olga on 22.01.17.
//

#include "SystemTools.h"

void SystemTools::showDotFile(string fileName) {
    string command = "dot -Tps " + fileName + "* -o tmp.ps";
    system(command.c_str());
    command = "gnome-open tmp.ps";
    system(command.c_str());
}