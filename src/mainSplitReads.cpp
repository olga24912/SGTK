#include <Logger/log_writers.hpp>
#include <ReadsSplitter/ReadsSplitter.h>
#include <ReadsSplitter/ReadsSplitter50.h>
#include <ReadsSplitter/SplitterByUnmappedEnd.h>

using namespace std;

int main(int argc, char **argv) {
    using namespace reads_splitter;
    logging::create_console_logger("../log.properties");

    if (argc != 5) {
        ERROR("expected 4 args: <0|1> - split50 or byUnmappedEnd, <RnafileName>, <ResFileName1>, <ResFileName2>");
        return 0;
    }

    if (std::string(argv[1]) == "0") {
        ReadsSplitter50 rs;
        rs.splitReads(argv[2], argv[3], argv[4]);
    } else if (std::string(argv[1]) == "1") {
        SplitterByUnmappedEnd rs;
        rs.splitReads(argv[2], argv[3], argv[4]);
    } else {
        ERROR("expected 0 or 1 on argv[1], but found " <<  argv[1]);
    }

    return 0;
}