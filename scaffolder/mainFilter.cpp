#include <Logger/log_writers.hpp>
#include "Filter/CommandParsers/Manager.h"

void create_console_logger(const std::string& log_props_file) {
    using namespace logging;

    //string log_props_file = cfg::get().log_filename;

    //if (!path::FileExists(log_props_file))
    //    log_props_file = path::append_path(dir, cfg::get().log_filename);

    logger *lg = create_logger(path::FileExists(log_props_file) ? log_props_file : "");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

int main(int argc, char **argv) {
    using namespace filter;
    create_console_logger("../log.properties");

    commands::Manager filter;
    filter.main();
    return 0;
}