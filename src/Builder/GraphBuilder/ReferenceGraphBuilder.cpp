#include "ReferenceGraphBuilder.h"

namespace builder {
    namespace graph_builder {
        void ReferenceGraphBuilder::evaluate() {
            INFO("Start build graph " + libName + " with type " + ContigGraph::Lib::typeToStr[libType]);
            if (tsvFileName == "") {
                alignmentREF(refContigFileName, queryContigsFileName);
            }

            generateVertex();

            if (tsvFileName == "") {
                createGraph(parseCoordFile("out.coords"));
            } else {
                createGraph(parseTSVFile(tsvFileName));
            }

            INFO("Finish build graph");
        }

        void ReferenceGraphBuilder::generateVertex() {
            INFO("start generateVertex");
            bool firstLib = (graph->getLibNum() == 1);

            seqan::SeqFileIn seqFileIn(queryContigsFileName.c_str());

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
            INFO("finish generateVertex");
        }

        void ReferenceGraphBuilder::createGraph(
                std::map<std::string, std::vector<ReferenceGraphBuilder::alignmentInfo>> contigsAlignment) {
            INFO("start createGraph");
            for (auto it = contigsAlignment.begin(); it != contigsAlignment.end(); ++it) {
                std::vector<alignmentInfo> contLPos = it->second;
                sort(contLPos.begin(), contLPos.end());
                for (int i = 0; i < (int) contLPos.size() - 1; ++i) {
                    int e = graph->incEdgeWeight(contigsId[contLPos[i].contigName],
                                                 contigsId[contLPos[i + 1].contigName],
                                                 contLPos[i].sr, contLPos[i].er, contLPos[i + 1].sr,
                                                 contLPos[i + 1].er);
                    graph->setEdgeChr(e, it->first);
                }
                for (int i = (int) contLPos.size() - 1; i > 0; --i) {
                    int e = graph->incEdgeWeight(contigsId[contLPos[i].contigName] ^ 1,
                                                 contigsId[contLPos[i - 1].contigName] ^ 1,
                                                 contLPos[i].er, contLPos[i].sr, contLPos[i - 1].er,
                                                 contLPos[i - 1].sr);
                    graph->setEdgeChr(e, it->first);
                }
            }
            INFO("finish createGraph");
        }

        std::map<std::string, std::vector<ReferenceGraphBuilder::alignmentInfo>>
        ReferenceGraphBuilder::parseCoordFile(std::string fileName) {
            INFO("start parse coord file=" << fileName);
            std::map<std::string, std::vector<ReferenceGraphBuilder::alignmentInfo>> contigsAlignment;
            std::ifstream in(fileName);

            int l, r, lq, rq, x;
            double xx;
            std::string rcont, qcont;

            while (in >> l >> r >> lq >> rq >> x >> x >> xx >> x >> x >> rcont >> qcont) {
                if (graph->getTargetLen(contigsId[qcont]) < minContigLen) {
                    continue;
                }

                if (lq > rq) {
                    qcont += "-rev";
                    std::swap(lq, rq);
                }

                if ((rq - lq) * 10 >= graph->getTargetLen(contigsId[qcont]) * 9) {
                    contigsAlignment[rcont].push_back(alignmentInfo(l, r, lq, rq, qcont));
                }
            }

            in.close();

            INFO("finish parse coord file");
            return contigsAlignment;
        }

        std::map<std::string, std::vector<ReferenceGraphBuilder::alignmentInfo>>
        ReferenceGraphBuilder::parseTSVFile(std::string fileName) {
            INFO("start parse TSV file=" << fileName);
            std::map<std::string, std::vector<ReferenceGraphBuilder::alignmentInfo>> contigsAlignment;
            std::ifstream in(fileName);

            std::string header;
            std::getline(in, header);

            int l, r, lq, rq, x;
            double xx;
            std::string bestGroup;

            std::string rcont, qcont;

            while (in >> l >> r >> lq >> rq >> rcont >> qcont >> x >> xx >> bestGroup) {
                std::string status;
                std::getline(in, status);
                std::getline(in, status);

                if (r - l < minContigLen) {
                    continue;
                }

                if (lq > rq) {
                    qcont += "-rev";
                    std::swap(lq, rq);
                }

                contigsAlignment[rcont].push_back(alignmentInfo(l, r, lq, rq, qcont));
            }

            in.close();

            INFO("finish parse TSV file");
            return contigsAlignment;
        }

        void ReferenceGraphBuilder::setMinContigLen(int minContigLen) {
            TRACE("setMinContigLen=" << minContigLen);
            this->minContigLen = minContigLen;
        }

        void ReferenceGraphBuilder::alignmentREF(std::string refFileName, std::string queryFileName) {
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
    }
}