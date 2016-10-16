#include <iostream>
#include "GraphControl.h"

using namespace std;

#include "GenTests/SplitReadsTest.h"

int main(int argc, char **argv) {
      SplitReadsTest splitReadsTest;

      splitReadsTest.genTest(string(argv[1]), "/home/olga/AU/bio-project/bio_scaffolder/test/ref.fasta",
                             "/home/olga/AU/bio-project/bio_scaffolder/test/reads.fasta", 100);

//    GraphControl graphControl;
//    graphControl.evaluate(argc, argv);
    return 0;
}