#include "ConnectionGraphBuilder.h"

namespace builder {
    namespace graph_builder {
        void ConnectionGraphBuilder::setConnectionFile(std::string file_name) {
            connectionFileName = file_name;
        }

        void ConnectionGraphBuilder::setContigFile(std::string file_name) {
            contigFileName = file_name;
        }

        void ConnectionGraphBuilder::evaluate() {
            initGraph();
            parseConnection();
        }

        void ConnectionGraphBuilder::initGraph() {
            bool firstLib = (graph->getLibNum() == 1);
            seqan::SeqFileIn seqFileIn(contigFileName.c_str());

            seqan::StringSet<seqan::CharString> ids;
            seqan::StringSet<seqan::Dna5String> seqs;

            seqan::readRecords(ids, seqs, seqFileIn);

            for (unsigned i = 0; i < length(ids); ++i) {
                std::string name = std::string(seqan::toCString(ids[i]));

                std::stringstream ss;
                ss << name;
                ss >> name;

                std::string seq = dna5ToString(seqan::toCString(seqs[i]), seqan::length(seqs[i]));

                contigsId[name] = 2 * i;
                contigsName.push_back(name);

                if (firstLib) {
                    graph->addVertex(2 * i, name, (int) seq.length());
                }

                name += "-rev";

                contigsId[name] = 2 * i + 1;
                contigsName.push_back(name);

                if (firstLib) {
                    graph->addVertex(2 * i + 1, name, (int) seq.length());
                }
            }
        }

        void ConnectionGraphBuilder::parseConnection() {
            std::ifstream in(connectionFileName);

            std::string fr, frd, sc, scd;
            double w;
            int len;

            while (in >> fr >> frd >> sc >> scd >> w >> len) {
                std::string info;
                getline(in, info);
                if (info != "") {
                    int i = 0;
                    while (i < info.size() && info[i] == ' ') {
                        ++i;
                    }
                    info = info.substr(i, info.size() - i);
                    while (info.size() > 0 && info[info.size() - 1] == ' ') {
                        info.resize(info.size() - 1);
                    }
                    if (info[0] == '"' && info[info.size() - 1] == '"') {
                        info = info.substr(1, info.size() - 2);
                    }
                }

                fr = fr.substr(1, fr.size() - 1);
                sc = sc.substr(1, sc.size() - 1);

                if (frd[0] == '-') {
                    fr += "-rev";
                }
                if (scd[0] == '-') {
                    sc += "-rev";
                }

                graph->addEdge(contigsId[fr], contigsId[sc], w, len, info);
                graph->addEdge(contigsId[sc]^1, contigsId[fr]^1, w, len, info);
            }

            in.close();
        }
    }
}
