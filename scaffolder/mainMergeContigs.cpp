#include <Logger/log_writers.hpp>
#include "ContigMerger/ContigMerger.h"

using namespace std;

void create_console_logger(const string& log_props_file) {
    using namespace logging;

    //string log_props_file = cfg::get().log_filename;

    //if (!path::FileExists(log_props_file))
    //    log_props_file = path::append_path(dir, cfg::get().log_filename);

    logger *lg = create_logger(path::FileExists(log_props_file) ? log_props_file : "");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

int main(int argc, char **argv) {
    create_console_logger("../log.properties");
    ContigMerger cm;
    cm.evaluate(string(argv[1]), string(argv[2]),
                string(argv[3]), "contigOUT.fasta", "samOUT.sam",
                string(argv[4]), string(argv[5]));

    return 0;
}