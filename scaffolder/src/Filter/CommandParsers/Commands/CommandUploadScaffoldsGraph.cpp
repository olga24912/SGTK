#include <seqan/file.h>
#include <seqan/store.h>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include "CommandUploadScaffoldsGraph.h"

namespace filter {
    namespace commands {
        void filter::commands::CommandUploadScaffoldsGraph::execute(std::string argv, filter::commands::State &state,
                                                                    filter::contig_graph::ContigGraph &graph) {

            using namespace contig_graph;
            std::stringstream ss(argv);
            std::string contigFileName;
            std::string scaffoldGraphFileName;

            ss >> contigFileName >> scaffoldGraphFileName;

            graph = ContigGraph();
            int libNum = graph.addLib("#ccffcc", "scaff", ContigGraph::Lib::SCAFF);

            addVertexis(graph, libNum, contigFileName);

            std::ifstream infoin(contigFileName);

            std::string cur;
            while (getline(infoin, cur)) {
                std::string u_name, v_name;
                int w, l;
                std::stringstream ss(cur);
                ss >> u_name >> v_name >> w >> l;

                int u_dir = (u_name[u_name.size() - 2] == '+');
                int v_dir = (v_name[v_name.size() - 2] == '+');

                u_name = u_name.substr(1, u_name.size() - 4);
                v_name = v_name.substr(1, v_name.size() - 4);

                int u = graph.getTargetId(u_name);
                int v = graph.getTargetId(v_name);
                if (u_dir == 0) u ^= 1;
                if (v_dir == 0) v ^= 1;

                graph.addEdge(u, v, libNum, w, 0, 0, 0, 0);
                graph.addEdge(v^1, u^1, libNum, w, 0, 0, 0, 0);
            }

            infoin.close();
        }

        void CommandUploadScaffoldsGraph::addVertexis(ContigGraph &graph, int libNum, std::string contigFileName) {
            seqan::SeqFileIn seqFileIn(contigFileName.c_str());
            seqan::StringSet<seqan::CharString> ids;
            seqan::StringSet<seqan::Dna5String> seqs;

            seqan::readRecords(ids, seqs, seqFileIn);

            for (unsigned i = 0; i < seqan::length(ids); ++i) {
                std::string contigName = getContigName(std::string(seqan::toCString(ids[i])));
                int seq_len = seqan::length(seqs[i]);

                graph.addVertex(i*2, contigName, seq_len);
                graph.addVertex(i*2 + 1, contigName + "-rev", seq_len);
            }
        }

        std::string CommandUploadScaffoldsGraph::getContigName(std::string s) {
            TRACE("getContigName s=" << s);
            std::string res = "";
            for (int i = 0; i < (int) s.size() && s[i] != ' '; ++i) {
                res += s[i];
            }
            return res;
        }
    }
}
