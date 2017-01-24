//
// Created by olga on 19.11.16.
//

#include "ReferenceGraphBuilder.h"

void ReferenceGraphBuilder::evaluate() {
    if (tsvFileName == "") {
        SystemAlignmentTools::alignmentREF(refContigFileName, queryContigsFileName);
    }

    generateVertex();

    if (tsvFileName == "") {
        createGraph(parseCoordFile("out.coords"));
    } else {
        createGraph(parseTSVFile(tsvFileName));
    }
}

void ReferenceGraphBuilder::generateVertex() {
    bool firstLib = (graph->getLibNum() == 1);

    SeqFileIn seqFileIn(queryContigsFileName.c_str());

    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;

    readRecords(ids, seqs, seqFileIn);

    for (unsigned i = 0; i < length(ids); ++i) {
        string name = string(toCString(ids[i]));

        stringstream ss;
        ss << name;
        ss >> name;

        cerr << name << endl;
        string seq = SeqanUtils::dna5ToString(toCString(seqs[i]), length(seqs[i]));

        contigsId[name] = 2 * i;
        contigsName.push_back(name);

        if (firstLib) {
            graph->addVertex(2*i, name, 0, seq.length());
        }

        name += "-rev";

        contigsId[name] = 2 * i + 1;
        contigsName.push_back(name);

        if (firstLib) {
            graph->addVertex(2 * i + 1, name, 0, seq.length());
        }
    }
}

string ReferenceGraphBuilder::getLibColor() {
    int color[3] = {255, 0, 0};
    return GraphUtils::colorToString(color);
}

void ReferenceGraphBuilder::createGraph(map<string, vector<ReferenceGraphBuilder::alignmentInfo>> contigsAlignment) {
    for (auto it = contigsAlignment.begin(); it != contigsAlignment.end(); ++it) {
        vector<alignmentInfo> contLPos = it->second;
        sort(contLPos.begin(), contLPos.end());
        for (int i = 0; i < (int)contLPos.size() - 1; ++i) {
            graph->incEdgeWeight(contigsId[contLPos[i].contigName], contigsId[contLPos[i + 1].contigName]);
        }
        for (int i = (int)contLPos.size() - 1; i > 0; --i) {
            graph->incEdgeWeight(contigsId[contLPos[i].contigName] ^ 1, contigsId[contLPos[i - 1].contigName] ^ 1);
        }
    }
}

map<string, vector<ReferenceGraphBuilder::alignmentInfo>> ReferenceGraphBuilder::parseCoordFile(string fileName) {
    map<string, vector<ReferenceGraphBuilder::alignmentInfo>> contigsAlignment;
    ifstream in(fileName);

    int l, r, lq, rq, x;
    double xx;
    string rcont, qcont;

    while (in >> l >> r >> lq >> rq >> x >> x >> xx >> x >> x >> rcont >> qcont) {
        if (graph->getTargetLength(contigsId[qcont]) < MIN_CONTIG) {
            continue;
        }

        if (lq > rq) {
            qcont += "-rev";
            swap(lq, rq);
        }

        if ((rq - lq) * 10 >= graph->getTargetLength(contigsId[qcont]) * 9) {
            contigsAlignment[rcont].push_back(alignmentInfo(l, r, lq, rq, qcont));
        }
    }

    in.close();

    return contigsAlignment;
}

map<string, vector<ReferenceGraphBuilder::alignmentInfo>> ReferenceGraphBuilder::parseTSVFile(string fileName) {
    cerr << "start parse TSV" << endl;
    map<string, vector<ReferenceGraphBuilder::alignmentInfo>> contigsAlignment;
    ifstream in(fileName);

    string header;
    getline(in, header);

    cerr << header << endl;
    int l, r, lq, rq, x;
    double xx;
    string bestGroup;

    string rcont, qcont;

    while (in >> l >> r >> lq >> rq >> rcont >> qcont >> x >> xx >> bestGroup) {
        string status;
        getline(in, status);
        getline(in, status);

        cerr << l << " " << r << " " << rcont << " " << qcont << endl;
        if (r - l < MIN_CONTIG) {
            continue;
        }

        if (lq > rq) {
            qcont += "-rev";
            swap(lq, rq);
        }

        contigsAlignment[rcont].push_back(alignmentInfo(l, r, lq, rq, qcont));
    }

    in.close();

    return contigsAlignment;
}