#ifndef SCAFFOLDER_GENEFINDER_H
#define SCAFFOLDER_GENEFINDER_H

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <Logger/logger.hpp>
#include <seqan/basic.h>
#include <seqan/gff_io.h>
#include "assert.h"


namespace filter {
    namespace alig_info {
        class GeneFinder {
        public:
            struct Gene {
                int b;
                int e;
                char strand;
                std::vector<Gene> exons;
                bool operator < (Gene sec) {
                    return b < sec.b;
                }
            };

        private:
            std::string gffFile;
            std::map<std::string, std::vector<Gene> > genesAnnotation;

        public:
            GeneFinder(std::string gff);
            Gene getGeneByCoord(std::string contigName, int coord);

        private:
            DECL_LOGGER("GeneFinder");
        };
    }
}


#endif //SCAFFOLDER_GENEFINDER_H
