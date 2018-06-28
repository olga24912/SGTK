#include <Logger/log_writers.hpp>
#include "Builder/GraphControl.h"

using namespace std;

int main(int argc, char **argv) {
    using namespace builder;
    logging::create_console_logger("/home/olga/bio-project/bio_scaffolder/scaffolder/src/log.properties");

    GraphControl graphControl;
    graphControl.evaluate(argc, argv);
    return 0;
}