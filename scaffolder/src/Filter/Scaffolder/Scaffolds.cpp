#include <ContigGraph/SeqanUtils.h>
#include "Scaffolds.h"
#include <seqan/seq_io.h>

void Scaffolds::print(std::string outFile) {
    std::SeqFileOut out(outFile.c_str());
    std::ofstream outInfo("out.info");

    std::StringSet<seqan::CharString> ids;
    std::StringSet<seqan::CharString> seqs;

    for (int i = 0; i < (int)scaffolds.size(); ++i) {
        if (scaffolds[i] != nullptr) {
            std::stringstream ss;
            ss << i;
            std::string readName = "path" + ss.str();
            outInfo << ">" + readName + " ";

            Node* cur = scaffolds[i];

            std::string seq = "";
            cur = cur->next;

            while (cur != nullptr) {
                if (cur->priv != nullptr) {
                    for (int j = 0; j < GAP_SIZE; ++j) {
                        seq += 'N';
                    }
                }

                if (((cur->id)&1) == 0) {
                    outInfo << "(" << contigsName[(cur->id)] << " "  << (cur->id)/2 << " +) ";
                } else {
                    outInfo << "(" << contigsName[(cur->id)^1] << " " << (cur->id)/2 << " -) ";
                }

                seq += contigs[cur->id];
                cur = cur->next;
            }

            outInfo << "\n";
            seqan::appendValue(ids, seqan::CharString(readName.c_str()));
            seqan::appendValue(seqs, seqan::Dna5String(seq));
        }
    }

    outInfo.close();
    seqan::writeRecords(out, ids, seqs);
}

void Scaffolds::saveContigs(std::string contigFile) {
    seqan::SeqFileIn seqFileIn(contigFile.c_str());
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;

    readRecords(ids, seqs, seqFileIn);

    contigsNode.resize(2 * seqan::length(ids));
    for (unsigned i = 0; i < seqan::length(ids); ++i) {
        std::string seq = SeqanUtils::dna5ToString(seqan::toCString(seqs[i]), seqan::length(seqs[i]));
        contigs.push_back(seq);
        contigs.push_back(createRevCompl(seq));

        contigsName.push_back(seqan::toCString(ids[i]));
        contigsName.push_back(seqan::toCString(ids[i]));

        contigsNode[2 * i].id = 2 * i;
        contigsNode[2 * i].rev = &contigsNode[2 * i + 1];
        contigsNode[2 * i + 1].id = 2 * i + 1;
        contigsNode[2 * i + 1].rev = &contigsNode[2 * i];
    }
}

std::string Scaffolds::createRevCompl(std::string s) {
    std::string res = s;
    reverse(res.begin(), res.end());
    for (int i = 0; i < (int)res.size(); ++i) {
        if (res[i] == 'A') {
            res[i] = 'T';
        } else if (res[i] == 'T') {
            res[i] = 'A';
        } else if (res[i] == 'C') {
            res[i] = 'G';
        } else if (res[i] == 'G') {
            res[i] = 'C';
        }
    }
    return res;
}

Scaffolds::Scaffolds(std::string contigFile) {
    saveContigs(contigFile);
}

void Scaffolds::addConnection(int id1, int id2) {
    Node* node1 = &contigsNode[id1];
    Node* node2 = &contigsNode[id2];

    node1->next = node2;
    scaffolds[id2] = nullptr;
}

void Scaffolds::brokeConnection(int id1) {
    Node* node1 = &contigsNode[id1];
    Node* node2 = node1->next;

    node1->next = nullptr;
    scaffolds[node2->id] = node2;
}