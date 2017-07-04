#include <stdio.h>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <sstream>
#include <Logger/log_writers.hpp>

void create_console_logger(const std::string& log_props_file) {
    using namespace logging;

    //string log_props_file = cfg::get().log_filename;

    //if (!path::FileExists(log_props_file))
    //    log_props_file = path::append_path(dir, cfg::get().log_filename);

    logger *lg = create_logger(path::FileExists(log_props_file) ? log_props_file : "");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

int main(int argc, char ** argv) {
    create_console_logger("../log.properties");
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
