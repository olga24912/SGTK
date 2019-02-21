#include "LongGraphBuilder.h"

void builder::graph_builder::LongGraphBuilder::setFileName(const std::string &fileName) {
    this->fileName = fileName;
}

void builder::graph_builder::LongGraphBuilder::initGraph() {
    seqan::CharString seqFileName(contigFileName.c_str());
    seqan::CharString id;
    seqan::Dna5String seq;

    seqan::SeqFileIn seqFileIn(toCString(seqFileName));
    int idd = 0;
    while (!seqan::atEnd(seqFileIn)) {
        readRecord(id, seq, seqFileIn);
        graph->addVertex(idd, std::string(seqan::toCString(id)), seqan::length(seq));
        idd += 1;
        graph->addVertex(idd, std::string(seqan::toCString(id)) + "-rev", seqan::length(seq));
        idd += 1;
    }
}

void builder::graph_builder::LongGraphBuilder::evaluate() {
    initGraph();
    std::ifstream in(fileName);
    std::string s;
    int i= 0;
    while (getline(in, s)) {
        updateEdges(s);
    }
    in.close();
}

void builder::graph_builder::LongGraphBuilder::updateEdges(std::string s) {
    std::stringstream ss(s);
    std::string readname;
    ss >> readname;
    if (aligInfo.size() > 0 && aligInfo[0].pacbioName != readname) {
        aligInfo.resize(0);
    }

    int len, cb, ce;
    char strand;
    std::string nodename;

    int nl, ncb, nce;
    ss >> len >> cb >> ce >> strand >> nodename >> nl >> ncb >> nce;

    if (!(ncb < 100 || nce > nl - 100 || (nce - ncb)*10 > 9*nl)) {
        return;
    }

    for (int i = 0; i < aligInfo.size(); ++i) {
        int u = -1, v = -1;
        if (aligInfo[i].cb < cb) {
            if ((aligInfo[i].ce - cb) * 10 < ce - aligInfo[i].cb) {
                u = graph->getTargetId(aligInfo[i].nodename);
                v = graph->getTargetId(nodename);
                if (aligInfo[i].strand == '-') {
                    u ^= 1;
                }
                if (strand == '-') {
                    v ^= 1;
                }
            }
        } else {
            if ((ce - aligInfo[i].cb) * 10 < aligInfo[i].ce - cb) {
                u = graph->getTargetId(nodename);
                v = graph->getTargetId(aligInfo[i].nodename);
                if (aligInfo[i].strand == '-') {
                    v ^= 1;
                }
                if (strand == '-') {
                    u ^= 1;
                }
            }
        }

        if (u != -1 && v != -1) {
            graph->incEdgeWeight(u, v, 0, 0, 0, 0);
            graph->incEdgeWeight(v ^ 1, u ^ 1, 0, 0, 0, 0);
        }
    }
    AligInfo alg(readname, len, cb, ce, strand, nodename);

    aligInfo.push_back(alg);
}

