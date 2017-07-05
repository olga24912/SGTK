#include <Logger/log_writers.hpp>
#include "Filter/CommandParsers/Manager.h"

int main(int argc, char **argv) {
    using namespace filter;
    logging::create_console_logger("../log.properties");

    commands::Manager filter;
    filter.main();
    return 0;
}