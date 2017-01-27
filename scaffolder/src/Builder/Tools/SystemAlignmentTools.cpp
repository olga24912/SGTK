//
// Created by olga on 10.10.16.
//

#include "SystemAlignmentTools.h"

void SystemAlignmentTools::alignmentRNA(string refFileName, string rnaFileName, string resFileName, std::string path) {
    std::string genomeDir = path + "/genomeDir";
    std::string command = "mkdir " + genomeDir;
    system(command.c_str());

    path += '/';

    command = "./STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ";
    command += genomeDir + " --genomeFastaFiles ";
    command += refFileName;
    command += " --genomeSAindexNbases 10";
    command += " --outFileNamePrefix " + path;
    //command += " --limitGenomeGenerateRAM 57888412374";

    cerr << command << endl;
    system(command.c_str());

    command = "./STAR --runThreadN 16 --genomeDir ";
    command += genomeDir + " --outReadsUnmapped Fastx --readFilesIn ";
    command += rnaFileName;
    command += " --outFileNamePrefix " + path;

    system(command.c_str());

    command = "mv " + path + "Aligned.out.sam ";
    command += path + resFileName;

    system(command.c_str());
}

void SystemAlignmentTools::alignmentREF(string refFileName, string queryFileName) {
    string command = "nucmer " + refFileName + " " + queryFileName;
    system(command.c_str());
    command = "show-coords out.delta -THrgl > out.coords";
    system(command.c_str());
}
