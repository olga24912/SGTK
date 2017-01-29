#include "SamFileWriteEdge.h"
#include <fstream>

SamFileWriteEdge::SamFileWriteEdge(std::string dir) {
    this->dir = dir + "/edges";

    std::string command = "mkdir " + this->dir;
    system(command.c_str());
}

void SamFileWriteEdge::writeEdge(int edgeID, seqan::BamAlignmentRecord read1, seqan::BamAlignmentRecord read2) {
    std::string fileName = getName(edgeID) + "#1.sam";

    std::ofstream out;
    out.open (fileName, std::ofstream::out | std::ofstream::app);
    seqan::BamFileOut fileOut(seqan::context(*fileIn), out, seqan::Sam());
    writeRecord(fileOut, read1);
    close(fileOut);
    out.close();


    fileName = getName(edgeID) + "#2.sam";
    out.open (fileName, std::ofstream::out | std::ofstream::app);
    seqan::BamFileOut fileOut2(seqan::context(*fileIn), out, seqan::Sam());
    writeRecord(fileOut2, read2);

    close(fileOut2);
    out.close();
}

std::string SamFileWriteEdge::getName(int edgeID) {
    std::string res = "";
    int x = edgeID;
    while (x != 0) {
        res += '0' + x%10;
        x /= 10;
    }
    std::reverse(res.begin(), res.end());
    res = dir + "/edge" + res;
    return res;
}
