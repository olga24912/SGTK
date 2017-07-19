#include <Filter/Statistics/WeightHistogram.h>
#include "CommandHistogram.h"

namespace filter {
    namespace commands {

        void CommandHistogram::execute(std::string argv, State &state, ContigGraph &graph) {
                std::stringstream ss(argv);

                int lib, step;
                ss >> lib >> step;

                INFO("Historgram lib=" << lib << " step=" << step);

                statistics::WeightHistogram::histogram(&graph, lib, step);
        }
    }
}