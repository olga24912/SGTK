//
// Created by olga on 31.10.16.
//

#include "SamFileWriteEdge.h"

SamFileWriteEdge::SamFileWriteEdge(string dir) {
    /*this->dir = dir;
    string commandRM = "rm -R ";
    string command = "mkdir ";
    commandRM.append(dir);
    command.append(dir);
    system(commandRM.c_str());
    system(command.c_str());*/
}

void SamFileWriteEdge::writeEdge(int edgeID, BamAlignmentRecord read1, BamAlignmentRecord read2) {
    /*string fileName = getName(edgeID) + "#1.sam";

    ofstream out;
    out.open (fileName, std::ofstream::out | std::ofstream::app);
    BamFileOut fileOut(context(*fileIn), out, Sam());
    writeRecord(fileOut, read1);
    close(fileOut);
    out.close();


    fileName = getName(edgeID) + "#2.sam";
    out.open (fileName, std::ofstream::out | std::ofstream::app);
    BamFileOut fileOut2(context(*fileIn), out, Sam());
    writeRecord(fileOut2, read2);

    close(fileOut2);
    out.close();*/
}

string SamFileWriteEdge::getName(int edgeID) {
    string res = "";
    int x = edgeID;
    while (x != 0) {
        res += '0' + x%10;
        x /= 10;
    }
    reverse(res.begin(), res.end());
    res = dir + "/edge" + res;
    return res;
}
