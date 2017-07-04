#include "SystemAlignmentTools.h"
#include <iostream>

void SystemAlignmentTools::alignmentRNA(std::string refFileName, std::string rnaFileName, std::string resFileName, std::string path) {
    INFO("start alignmentRNA refFileName=" << refFileName << " rnaFileName=" << rnaFileName
                                           << " resFileName=" << resFileName << "outputPath=" << path);

    std::string genomeDir = path + "/genomeDir";
    std::string command = "mkdir " + genomeDir;
    system(command.c_str());
    DEBUG("command: " << command);

    path += '/';

    command = "./STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ";
    command += genomeDir + " --genomeFastaFiles ";
    command += refFileName;
    command += " --genomeSAindexNbases 10";
    command += " --outFileNamePrefix " + path;
    command += " --limitGenomeGenerateRAM 90000000000";
    system(command.c_str());
    DEBUG("command: " << command);

    command = "./STAR --runThreadN 16 --genomeDir ";
    command += genomeDir + " --outReadsUnmapped Fastx --readFilesIn ";
    command += rnaFileName;
    command += " --outFileNamePrefix " + path;
    system(command.c_str());
    DEBUG("command: " << command);

    command = "mv " + path + "Aligned.out.sam ";
    command += path + resFileName;
    system(command.c_str());
    DEBUG("command: " << command);
    INFO("finish alignmentRNA");
}

void SystemAlignmentTools::alignmentREF(std::string refFileName, std::string queryFileName) {
    INFO("start alignmentREF refFileName=" << refFileName << " queryFileName=" << queryFileName);

    std::string command = "nucmer " + refFileName + " " + queryFileName;
    system(command.c_str());
    DEBUG("command: " << command);

    command = "show-coords out.delta -THrgl > out.coords";
    std::cerr << command << std::endl;
    system(command.c_str());
    DEBUG("command: " << command);

    INFO("finish alignmentREF");
}
