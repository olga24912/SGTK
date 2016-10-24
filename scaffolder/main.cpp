#include <iostream>
#include "GraphControl.h"

using namespace std;

#include "GenTests/SplitReadsTest.h"
#include "ContigGraph/InteractiveFilter.h"

int main(int argc, char **argv) {
      /*SplitReadsTest splitReadsTest;

      splitReadsTest.genTest(string(argv[1]), "/home/olga/AU/bio-project/bio_scaffolder/test/ref.fasta",
                             "/home/olga/AU/bio-project/bio_scaffolder/test/reads.fasta", 100, 11);
*/
    /*GraphControl graphControl;
    graphControl.evaluate(argc, argv);*/
    InteractiveFilter::main();
    return 0;
}