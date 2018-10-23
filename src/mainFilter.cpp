/*
 * Graph filtering main
 */

#include <Logger/log_writers.hpp>
#include "Filter/CommandParsers/Manager.h"

int main(int argc, char **argv) {
    using namespace filter;
    logging::create_console_logger("");

    commands::Manager filter;
    if (argc > 1) {
        filter.main(std::string(argv[1]));
    } else {
        filter.main("");
    }
    return 0;
}