#ifndef SCAFFOLDER_FILEVALIDATOR_H
#define SCAFFOLDER_FILEVALIDATOR_H

#include <vector>
#include <Filter/Filters/Filter.h>

namespace filter {
    namespace writers {
        class FileValidator {
        public:
            virtual bool isGoodVertexSet(std::vector<int> vert, Filter *filter) {
                return true;
            }

        protected:
            DECL_LOGGER("FileValidator");
        };
    }
}

#endif //SCAFFOLDER_FILEVALIDATOR_H
