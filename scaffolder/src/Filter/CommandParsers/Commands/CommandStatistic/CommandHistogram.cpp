#include <Filter/Statistics/WeightHistogram.h>
#include "CommandHistogram.h"

namespace filter {
    namespace commands {

        void CommandHistogram::execute(std::string argv, State &state, ContigGraph *filter) {
                std::stringstream ss(argv);

                int lib, step;
                ss >> lib >> step;

                INFO("Historgram lib=" << lib << " step=" << step);

                statistics::WeightHistogram::histogram(filter, lib, step);
        }
    }
}