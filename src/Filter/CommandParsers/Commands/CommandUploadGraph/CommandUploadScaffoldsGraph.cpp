#include <seqan/file.h>
#include <seqan/store.h>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <Filter/ContigGraph/ContigGraph.h>
#include "CommandUploadScaffoldsGraph.h"

namespace filter {
    namespace commands {
        void filter::commands::CommandUploadScaffoldsGraph::execute(std::string argv, filter::commands::State &state,
                                                                    filter::contig_graph::ContigGraph &graph) {

            using namespace contig_graph;
            INFO("start upload scaffold graph");
            std::stringstream ss(argv);
            std::string contigFileName;
            std::string scaffoldGraphFileName;

            ss >> contigFileName >> scaffoldGraphFileName;

            graph = ContigGraph();
            int libNum = graph.addLib("#00ffff", "scaff", ContigGraph::Lib::SCAFF);

            addVertexis(graph, libNum, contigFileName);

            std::ifstream infoin(scaffoldGraphFileName);

            std::string cur;
            while (getline(infoin, cur)) {
                std::string u_name, v_name;
                int w, l;
                std::stringstream ss(cur);
                ss >> u_name >> v_name >> w >> l;

                std::cerr << u_name << " " << v_name << " " << w << " " << l << "\n";

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
            INFO("finish upload scaffold graph");
        }

        void CommandUploadScaffoldsGraph::addVertexis(ContigGraph &graph, int libNum, std::string contigFileName) {
            DEBUG("start add vert");
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
            DEBUG("finish add vert");
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
