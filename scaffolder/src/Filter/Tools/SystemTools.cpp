#include "SystemTools.h"

namespace filter {
    namespace tools {
        void SystemTools::showDotFile(std::string fileName) {
            std::string command = "dot -Tps " + fileName + "* -o tmp.ps";
            system(command.c_str());
            command = "gnome-open tmp.ps";
            system(command.c_str());
        }
    }
}