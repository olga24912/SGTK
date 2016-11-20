//
// Created by olga on 19.11.16.
//

#include "ReferenceGraphBuilder.h"

void ReferenceGraphBuilder::evaluate() {
    WorkWithOtherTools::alignmentREF(refContigFileName, queryContigsFileName);

    generateVertex();
    createGraph("out.coords");
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

void ReferenceGraphBuilder::createGraph(string fileName) {
    ifstream in(fileName);

    int l, r, lq, rq, x;
    double xx;
    string rcont, qcont;

    string currentRef = "";
    vector< pair<int, string> > contLPos;
    while (in >> l >> r >> lq >> rq >> x >> x >> xx >> x >> x >> rcont >> qcont) {
        if (graph->getTargetLength(contigsId[qcont]) < MIN_CONTIG) {
            continue;
        }
        if (rcont != currentRef) {
            sort(contLPos.begin(), contLPos.end());
            for (int i = 0; i < (int)contLPos.size() - 1; ++i) {
                graph->incEdgeWeight(contigsId[contLPos[i].second], contigsId[contLPos[i + 1].second]);
            }
            for (int i = (int)contLPos.size() - 1; i > 0; --i) {
                graph->incEdgeWeight(contigsId[contLPos[i].second] ^ 1, contigsId[contLPos[i - 1].second] ^ 1);
            }

            contLPos.resize(0);
            currentRef = rcont;
        }
        if (lq > rq) {
            qcont += "-rev";
            swap(lq, rq);
        }

        cerr << rq - lq << " " << graph->getTargetLength(contigsId[qcont]) << " " << contigsId[qcont] << " " << qcont << endl;
        if ((rq - lq) * 10 >= graph->getTargetLength(contigsId[qcont]) * 9) {
            cerr << "push" << endl;
            contLPos.push_back(make_pair(l, qcont));
        }
    }

    sort(contLPos.begin(), contLPos.end());
    for (int i = 0; i < (int)contLPos.size() - 1; ++i) {
        graph->incEdgeWeight(contigsId[contLPos[i].second], contigsId[contLPos[i + 1].second]);
    }

    for (int i = (int)contLPos.size() - 1; i > 0; --i) {
        graph->incEdgeWeight(contigsId[contLPos[i].second] ^ 1, contigsId[contLPos[i - 1].second] ^ 1);
    }

    in.close();
}
