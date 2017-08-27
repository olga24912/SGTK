#ifndef SCAFFOLDER_RULECOVERING_H
#define SCAFFOLDER_RULECOVERING_H

#include <Filter/CommandParsers/Commands/Command.h>
#include "Rule.h"
#include <seqan/file.h>
#include <seqan/bam_io.h>

namespace filter {
    namespace scaffolder {
        using namespace commands;
        using namespace statistics;
        class RuleCovering : public Rule {
        private:
            std::vector<State::BamFiles> bamFiles;
            InfoAboutContigsAlig* alig;

            void checkThisEdge(ContigGraph *graph, int e);
            double getCover1(ContigGraph *graph, State::BamFiles bamFile, int e);
            double getCover2(ContigGraph *graph, State::BamFiles bamFile, int e);
        public:
            void simplifyGraph(ContigGraph *filter) override;

            void setBamFiles(std::vector<State::BamFiles> bamFiles) {
                this->bamFiles = bamFiles;
            }

            void setAligFile(InfoAboutContigsAlig* alig) {
                this->alig = alig;
            }

            int getCover(int cb, int ce, seqan::BamFileIn &inFile, const seqan::BamIndex<seqan::Bai> &baiIndex, int rID,
                 int isRev) const;
        };
    }
}


#endif //SCAFFOLDER_RULECOVERING_H
