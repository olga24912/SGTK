/*
 * Graph building main
 */

#include <Logger/log_writers.hpp>
#include "Builder/GraphControl.h"

int main(int argc, char **argv) {
    using namespace builder;
    logging::create_console_logger("");

    GraphControl graphControl;
    graphControl.evaluate(argc, argv);
    return 0;
}