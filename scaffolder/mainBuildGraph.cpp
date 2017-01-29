#include "Builder/GraphControl.h"

using namespace std;

int main(int argc, char **argv) {
    GraphControl graphControl;
    graphControl.evaluate(argc, argv);

   /* ContigMerger cm;
    cm.evaluate(string(argv[1]), string(argv[2]),
    string(argv[3]), "contigOUT.fasta", "samOUT.sam",
    string(argv[4]), string(argv[5]));*/
    return 0;
}