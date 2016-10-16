//
// Created by olga on 08.10.16.
//

#include "ReadsAlignment.h"
#include "Tools/FastaToolsOut.h"
#include "Tools/FastaToolsIn.h"

void ReadsAlignment::alignment(string refFileName, unordered_map<string, string> fullReads,
                               string alignfileName, string resFile1, string resFile2) {
    FastaToolsOut ftout1;
    FastaToolsOut ftout2;

    ftout1.putFileName(resFile1);
    ftout2.putFileName(resFile2);


    unordered_map<string, string> contigs = getContigs(refFileName);
    unordered_map<string, string> pair_contigs; //TODO
    vector<string> contigName;

    BamFileIn bamFile;
    open(bamFile, alignfileName.c_str());

    BamHeader samHeader;
    readHeader(samHeader, bamFile);

    typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;

    TBamContext const &bamContext = context(bamFile);
    size_t contig_num = length(contigNames(bamContext));

    string name;
    for (int i = 0; i < static_cast<int>(contig_num); ++i) {
        name = string(toCString(contigNames(bamContext)[i]));
        contigName.push_back(name);
    }

    BamAlignmentRecord read;
    while (!atEnd(bamFile)) {
        readRecord(read, bamFile);
        string name = string(toCString(read.qName));
        cerr << "? " << name << endl;
        if (read.rID != -1) {
            string cutName = SeqanUtils::cutReadName(read);

            int cnt = -1;
            cerr << cutName << endl;
            if (name[name.size() - 1] == '1' && !hasFlagRC(read)) {
                cerr << read.rID << endl;
                cerr << getAlignmentLengthInRef(read) << endl;
                cnt = alignFull(contigs[contigName[read.rID]], fullReads[cutName], read.beginPos);
            } else if (name[name.size() - 1] == '2' && !hasFlagRC(read)) {
                string ref = contigs[contigName[read.rID]];
                int bg = (int)ref.size() - (read.beginPos + getAlignmentLengthInRef(read));
                reverse(ref.begin(), ref.end());
                string rd = fullReads[cutName];
                reverse(rd.begin(), rd.end());

                cerr << read.rID << endl;
                cerr << getAlignmentLengthInRef(read) << endl;
                cerr << rd << endl;
                cerr << bg << endl;

                cnt = (int)rd.size() - alignFull(ref, rd, bg);
            }

            cerr << cnt << endl;
            string rseq = fullReads[cutName];
            cerr << rseq.size() << endl;
            if (cnt != -1 && cnt != 0 && cnt != rseq.size()) {
                ftout1.write(cutName.append("/1"), rseq.substr(0, cnt));
                cutName.resize(cutName.size() - 2);
                ftout2.write(cutName.append("/2"), rseq.substr(cnt, fullReads[cutName].size() - cnt));
            }
        }
    }
    close(bamFile);
    ftout1.close();
    ftout2.close();
}

unordered_map<string, string> ReadsAlignment::getContigs(string fileName) {
    unordered_map<string, string> contig;
    FastaToolsIn ftin;
    ftin.parse(fileName);

    while (ftin.next()) {
        contig[ftin.currentName()] = ftin.currentRef();
    }

    ftin.close();
    return contig;
}

int ReadsAlignment::alignFull(string ref, string read, int pos) {
    ref = ref.substr(pos, ref.size() - pos);
    ref.resize(min(ref.size(), read.size() * 3));

    vector< vector<int> > d(ref.size() + 1, vector<int> (read.size() + 1, -inf));
    d[0][0] = 0;
    for (int i = 1; i <= read.size(); ++i) {
        d[0][i] = -i;
    }

    for (int i = 1; i <= ref.size(); ++i) {
        d[i][0] = -i;
    }

    int x = 0, y = 0;

    for (int i = 1; i <= ref.size(); ++i) {
        for (int j = 1; j <= read.size(); ++j) {
            if (ref[i - 1] == read[j - 1]) {
                d[i][j] = max(d[i][j], d[i - 1][j - 1] + 1);
            }
            d[i][j] = max(d[i][j], d[i - 1][j - 1] - 1);
            d[i][j] = max(d[i][j], d[i - 1][j] - 1);
            d[i][j] = max(d[i][j], d[i][j - 1] - 1);
            if (d[i][j] > d[x][y] || (d[i][j] == d[x][y] && j < y)) {
                x = i;
                y = j;
            }
        }
    }

    return y;
}
