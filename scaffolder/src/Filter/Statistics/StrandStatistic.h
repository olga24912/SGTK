#ifndef SCAFFOLDER_STRANDSTATISTIC_H
#define SCAFFOLDER_STRANDSTATISTIC_H


#include <Filter/ContigGraph/ContigGraph.h>
#include <Filter/AligInfo/GeneAnnotationAns.h>
#include <Filter/AligInfo/GeneAnnotationContigs.h>
#include "Statistic.h"

namespace filter {
    namespace statistics {
        using namespace contig_graph;

        class StrandStatistic : public Statistic {
        private:
            std::string rightout = "statistic/rightStrandRightLock";
            std::string rightstrandWrongLock = "statistic/rightStrandWrongLock";
            std::string wrongstrandRightLock = "statistic/wrongStrandRightLock";
            std::string wrongout = "statistic/wrongStrandWrongLock";
            std::string noneout = "statistic/geneNotFound";

            std::ofstream outrht;
            std::ofstream outrhtwng;
            std::ofstream outwngrht;
            std::ofstream outwng;
            std::ofstream outnone;

            alig_info::GeneAnnotationAns geneans;
            alig_info::GeneAnnotationContigs geneContig;
            alig_info::InfoAboutContigsAlig aligInfo;

            filter::contig_graph::ContigGraph *graph;
        public:

            StrandStatistic(filter::contig_graph::ContigGraph *graph, alig_info::InfoAboutContigsAlig aligInfo,
                            std::string gffFileRef, std::string gffFileContig, std::string out): aligInfo(aligInfo),
                    graph(graph), geneans(gffFileRef, aligInfo), geneContig(gffFileContig) {
                outrht.open(rightout);
                outrhtwng.open(rightstrandWrongLock);
                outwngrht.open(wrongstrandRightLock);
                outwng.open(wrongout);
                outnone.open(noneout);
            }

            void calculateStatistic();

            void calculateStatisticForOneEdge(int e);
            void calculateStatisticForOneEdge1(int e);
            void calculateStatisticForOneEdge2(int e);

            void printGene(alig_info::GeneFinder::Gene gene, std::ofstream &outs);

            bool printAligInfo(alig_info::InfoAboutContigsAlig alig, ContigGraph *graph, int e, std::ofstream &outs);

            void printInfo(filter::contig_graph::ContigGraph *graph,
                                                                const filter::alig_info::InfoAboutContigsAlig &aligInfo, std::ofstream &outs, int v,
                                                                int e, const filter::alig_info::GeneFinder::Gene &gener,
                                                                const filter::alig_info::GeneFinder::Gene &genec);

            ~StrandStatistic() {
                outrht.close();
                outrhtwng.close();
                outwngrht.close();
                outwng.close();
                outnone.close();
            }
        };
    }
}



#endif //SCAFFOLDER_STRANDSTATISTIC_H
