//
// Created by olga on 10.10.16.
//

#include "WorkWithOtherTools.h"

void WorkWithOtherTools::alignmentRNA(string refFileName, string rnaFileName, string resFileName) {
    system("mkdir genomeDir");

    string command = "./STAR --runThreadN 16 --runMode genomeGenerate --genomeDir genomeDir --genomeFastaFiles ";
    command += refFileName;
    command += " --genomeSAindexNbases 10";
    //command += " --limitGenomeGenerateRAM 57888412374";

    cerr << command << endl;
    system(command.c_str());

    command = "./STAR --runThreadN 16 --genomeDir genomeDir --outReadsUnmapped Fastx --readFilesIn ";
    command += rnaFileName;

    system(command.c_str());

    command = "mv Aligned.out.sam ";
    command += resFileName;

    system(command.c_str());
}

void WorkWithOtherTools::alignmentREF(string refFileName, string queryFileName) {
    string command = "nucmer " + refFileName + " " + queryFileName;
    system(command.c_str());
    command = "show-coords out.delta -THrgl > out.coords";
    system(command.c_str());
}
