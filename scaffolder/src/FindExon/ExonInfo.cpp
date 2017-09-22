#include <seqan/bam_io.h>
#include <Logger/logger.hpp>
#include "ExonInfo.h"

namespace findExon {

    void findExon::ExonInfo::writeExonBlock(std::ostream &out) {
        out << contigName << " ";

        for (int i = 0; i < exons.size(); ++i) {
            out << "( " << exons[i].b << " " << exons[i].e << " " << exons[i].cover << " " << exons[i].genId << " ) ";
        }
        out << "\n";
    }

    void findExon::ExonInfo::addReadInfo(seqan::BamAlignmentRecord read) {
        int bgPos = read.beginPos;

        seqan::String<seqan::CigarElement<> > cigar = read.cigar;
        for (seqan::CigarElement<> elem : cigar) {
            if (elem.operation == 'M' || elem.operation == 'D' || elem.operation == '=' || elem.operation == 'X') {
                cover[bgPos] += 1;
                cover[bgPos + elem.count] -= 1;
                bgPos += elem.count;
            } else if (elem.operation == 'N') {
                misCover[bgPos] += 1;
                misCover[bgPos + elem.count] -= 1;
                bgPos += elem.count;
            }
        }
    }

    void ExonInfo::finish() {
        DEBUG("finish contig " << contigName);
        const int MAX_SPACE_SIZE = 100;
        const int MIN_COV = 2;
        const int MIN_MIS_COV = 3;

        int cur = 0;
        int pos = 0;
        Exon exon;
        exon.b = -1;
        exon.e = -1;
        exon.genId = 0;
        exon.cover = 0;
        int sum = 0;

        std::vector<int> next(len + 1, -1);
        for (int i = 0; i <= len; ++i) {
            cur += cover[i];
            cover[i] = cur;
        }

        for (int i = len; i >= 0; --i) {
            if (cover[i] >= MIN_COV) {
                next[i] = i;
            } else if (i < len) {
                next[i] = next[i + 1];
            }
        }

        TRACE("calcul next");

        while (pos <= len) {
            if (cover[pos] >= MIN_COV) {
                if (exon.b == -1) {
                    exon.b = pos;
                }
                exon.e = pos;
                sum += cover[pos];
            } else {
                if (next[pos] != -1 && next[pos] - pos <= MAX_SPACE_SIZE && exon.b != -1) {
                    exon.e = pos;
                    sum += cover[pos];
                } else if (exon.b != -1) {
                    exon.cover = sum*1./(exon.e - exon.b);
                    if (exon.e - exon.b >= MAX_SPACE_SIZE && exon.cover > 10) {
                        exons.push_back(exon);
                    }
                    exon.b = -1;
                    exon.e = -1;
                    exon.cover = 0;
                    sum = 0;
                }
            }
            cover[pos] = 0;
            ++pos;
        }

        if (pos == len + 1 && exon.b != -1) {
            exon.cover = sum*1./(exon.e - exon.b);
            exons.push_back(exon);
        }

        findExonCon();
    }

    void ExonInfo::findExonCon() {
        TRACE("start find exon con");
        const int MIN_MIS_COV = 2;
        int cur = 0;

        std::vector<int> next(len + 1, -1);
        for (int i = 0; i <= len; ++i) {
            cur += misCover[i];
            misCover[i] = cur;
        }

        for (int i = len; i >= 0; --i) {
            if (misCover[i] >= MIN_MIS_COV) {
                next[i] = i;
            } else if (i < len) {
                next[i] = next[i + 1];
            }
            misCover[i] = 0;
        }

        int id = 0;
        for (int i = 1; i < exons.size(); ++i) {
            if (hasConnection(exons[i - 1], exons[i], next)) {
                exons[i].genId = id;
            } else {
                ++id;
                exons[i].genId = id;
            }
        }
        TRACE("finish find exon con");
    }

    bool ExonInfo::hasConnection(ExonInfo::Exon &exon, ExonInfo::Exon &exon1, std::vector<int> &next) {
        const int MAX_SPACE_SIZE = 100;
        int pos = exon.e;
        while (pos + MAX_SPACE_SIZE < exon1.b) {
            if (next[pos] == -1) {
                return false;
            }
            if (next[pos] - pos > MAX_SPACE_SIZE) {
                return false;
            }
            if (pos == next[pos]) pos += 1;
            else pos = next[pos];

        }

        return true;
    }
}