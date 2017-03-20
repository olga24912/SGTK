#include "SystemAlignmentTools.h"
#include <iostream>

void SystemAlignmentTools::alignmentRNA(std::string refFileName, std::string rnaFileName, std::string resFileName, std::string path) {
    std::string genomeDir = path + "/genomeDir";
    std::string command = "mkdir " + genomeDir;
    system(command.c_str());

    path += '/';

    command = "./STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ";
    command += genomeDir + " --genomeFastaFiles ";
    command += refFileName;
    command += " --genomeSAindexNbases 10";
    command += " --outFileNamePrefix " + path;
    command += " --limitGenomeGenerateRAM 90000000000";

    std::cerr << command << std::endl;
    system(command.c_str());

    command = "./STAR --runThreadN 16 --genomeDir ";
    command += genomeDir + " --outReadsUnmapped Fastx --readFilesIn ";
    command += rnaFileName;
    command += " --outFileNamePrefix " + path;

    std::cerr << command << std::endl;
    system(command.c_str());

    command = "mv " + path + "Aligned.out.sam ";
    command += path + resFileName;

    system(command.c_str());
}

void SystemAlignmentTools::alignmentREF(std::string refFileName, std::string queryFileName) {
    std::string command = "nucmer " + refFileName + " " + queryFileName;
    std::cerr << command << std::endl;
    system(command.c_str());
    command = "show-coords out.delta -THrgl > out.coords";
    std::cerr << command << std::endl;
    system(command.c_str());
}
