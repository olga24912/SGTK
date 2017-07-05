#include <stdio.h>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <sstream>
#include <Logger/log_writers.hpp>

int main(int argc, char ** argv) {
    logging::create_console_logger("../log.properties");
    FILE * in = fopen(argv[1], "r");
    std::string pr = std::string(argv[2]);
    int threadN = atoi(argv[3]);
    int lineCnt = atoi(argv[4]);

    int mc = lineCnt/threadN + 1;

    char str[150];
    str[0] = 0;
    fgets(str, 150, in);
    char fr = str[0];
    for (int i = 0; i < threadN; ++i) {
        std::stringstream ss;
        ss << pr << i;
        std::string fileName(ss.str());

        FILE* file = fopen(fileName.c_str(), "w");

        int cur = 1;
        int oc = cur;
        while (cur < mc) {
            fputs(str, file);
            while (fgets(str, 150, in) != NULL && str[0] != fr) {
                fputs(str, file);
                ++cur;
            }
            if (oc == cur) break;
            oc = cur;
        }

        fclose(file);
    }

    fclose(in);

    return 0;
}
