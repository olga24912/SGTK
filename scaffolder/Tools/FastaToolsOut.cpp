//
// Created by olga on 16.10.16.
//

#include "FastaToolsOut.h"

void FastaToolsOut::close() {
    fout.close();
}

void FastaToolsOut::putFileName(string fileName) {
    this->fileName = fileName;

    fout.open(fileName);
}

void FastaToolsOut::write(string name, string ref) {
    fout << ">" << name << "\n";
    fout << ref<< "\n";
}
