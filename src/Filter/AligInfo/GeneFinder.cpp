#include "GeneFinder.h"

namespace filter {
    namespace alig_info {
        GeneFinder::GeneFinder(std::string gff) {
            seqan::GffFileIn gffIn(gff.c_str());

            seqan::GffRecord record;
            while (!seqan::atEnd(gffIn)) {
                try {
                    seqan::readRecord(record, gffIn);
                    int b = record.beginPos, e = record.endPos;
                    std::string contigName = std::string(seqan::toCString(record.ref));

                    Gene ex;
                    ex.b = b;
                    ex.e = e;
                    ex.strand = record.strand;

                    if (std::string(seqan::toCString(record.type)) == "gene") {
                        genesAnnotation[contigName].push_back(ex);
                    } else if (std::string(seqan::toCString(record.type)) == "CDS") {
                        genesAnnotation[contigName][genesAnnotation[contigName].size() - 1].exons.push_back(ex);
                    }
                } catch (seqan::Exception const & e) {
                    WARN(e.what());
                }
            }
        }

        GeneFinder::Gene GeneFinder::getGeneByCoord(std::string contigName, int coord) {
            DEBUG("get gene by coord " << contigName << " " << coord);
            std::vector<Gene> genes = genesAnnotation[contigName];

            Gene fg;
            fg.e = coord;
            fg.b = coord;

            int pos = std::lower_bound(genes.begin(), genes.end(), fg) - genes.begin();

            if (pos > 0) {
                --pos;
                if (genes[pos].b  - 5 <= coord && genes[pos].e + 5 >= coord) {
                    return genes[pos];
                }
                ++pos;
            }

            if (pos + 1 < genes.size()) {
                ++pos;
                if (genes[pos].b  - 5 <= coord && genes[pos].e + 5 >= coord) {
                    return genes[pos];
                }
                --pos;
            }

            if (pos == genes.size()) {
                Gene ret;
                ret.b = -1;
                return ret;
            } else {
                return genes[pos];
            }
        }
    }
}