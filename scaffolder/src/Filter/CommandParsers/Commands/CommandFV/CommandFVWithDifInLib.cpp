#include <Filter/Writers/FileValidator/ValidatorWithDifInLib.h>
#include "CommandFVWithDifInLib.h"

namespace filter {
    namespace commands {

        void CommandFVWithDifInLib::setFV(State &state, std::string argv) {
            INFO("set file validator with dif in libs " << argv);

            int libNum;
            std::vector<int> libs;
            std::stringstream ss(argv);
            while (ss >> libNum) {
                libs.push_back(libNum);
            }

            state.validator = new writers::ValidatorWithDifInLib(libs);
        }
    }
}
