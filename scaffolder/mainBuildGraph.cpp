#include "Builder/GraphControl.h"

using namespace std;


int main(int argc, char **argv) {
      /*SplitReadsTest splitReadsTest;

      splitReadsTest.genTest(string(argv[1]), "/home/olga/AU/bio-project/bio_scaffolder/test/ref.fasta",
                             "/home/olga/AU/bio-project/bio_scaffolder/test/reads.fasta", 100, 11);
*/
    GraphControl graphControl;
    graphControl.evaluate(argc, argv);

    //InteractiveFilter::main();

   /* ContigMerger cm;
    cm.evaluate(string(argv[1]), string(argv[2]),
    string(argv[3]), "contigOUT.fasta", "samOUT.sam",
    string(argv[4]), string(argv[5]));*/
    return 0;
}