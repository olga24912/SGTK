#ifndef SCAFFOLDER_GRAPHCONTROL_H
#define SCAFFOLDER_GRAPHCONTROL_H

#include "ContigGraph/ContigGraph.h"
#include "Builder/GraphBuilder/DNAPairReadGraphBuilder.h"
#include "Builder/ExternalLib/optionparser.h"
#include "Builder/GraphBuilder/PairReadGraphBuilder.h"
#include "Builder/GraphBuilder/RNASplitReadGraphBuilder.h"

// that class parse args and call other class for
// build conection and keep info about conection.
class GraphControl {
private:
    ContigGraph graph;
    std::string path = ".";

    struct Arg: public option::Arg
    {
        static void printError(const char* msg1, const option::Option& opt, const char* msg2)
        {
            fprintf(stderr, "%s", msg1);
            fwrite(opt.name, opt.namelen, 1, stderr);
            fprintf(stderr, "%s", msg2);
        }

        static option::ArgStatus Unknown(const option::Option& option, bool msg)
        {
            if (msg) printError("Unknown option '", option, "'\n");
            return option::ARG_ILLEGAL;
        }

        static option::ArgStatus Required(const option::Option& option, bool msg)
        {
            if (option.arg != 0)
                return option::ARG_OK;

            if (msg) printError("Option '", option, "' requires an argument\n");
            return option::ARG_ILLEGAL;
        }

        static option::ArgStatus Numeric(const option::Option& option, bool msg)
        {
            char* endptr = 0;
            if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
            if (endptr != option.arg && *endptr == 0)
                return option::ARG_OK;

            if (msg) printError("Option '", option, "' requires a numeric argument\n");
            return option::ARG_ILLEGAL;
        }


        static option::ArgStatus NewLib(const option::Option& option, bool msg)
        {
            std::string arg = option.arg;
            if (arg == "RNA_PAIR" || arg == "RNA_SPLIT" || arg == "DNA_PAIR" || arg == "REF") {
                return option::ARG_OK;
            }

            return option::ARG_ILLEGAL;
        }
    };

    enum  optionIndex { UNKNOWN, HELP, NEW, SAM1, SAM2, MINCONTIGLEN,
        DIST_BETWEEN_PAIR_READS, REFFILE, QUERYFILE, READSFILE, LIBNAME, TSVFILE, OUTPUT};
    const option::Descriptor usage[14] = {
            { UNKNOWN, 0,"", "",        Arg::Unknown, "USAGE: scaffolder [options]\n\n"
                                                              "Options:" },
            { HELP,         0,"", "help",         Arg::None,    "  \t--help  \tPrint usage and exit." },
            { NEW,          0,"n","new",          Arg::NewLib,  "  -n <RNA_PAIR|RNA_SPLIT|DNA_PAIR|REF>, \t--new=<RNA_PAIR|RNA_SPLIT|DNA_PAIR|REF>"
                                                              "  \tNeed to take one of string RNA_PAIR, RNA_SPLIT, DNA_PAIR" },
            { SAM1,         0,"f","sam1",         Arg::Required,"  -f <arg>, \t--sam1=<arg>  \tMust have an argument,"
                                                              "file name of first sam file. def = read1.sam" },
            { SAM2,         0,"s","sam2",         Arg::Required,"  -s <arg>, \t--sam2=<arg>  \tMust have an argument,"
                                                              "file name of second sam file. def = read2.sam" },
            { MINCONTIGLEN, 0,"c","mincontiglen", Arg::Numeric, "  -c <num>, \t--mincontiglen=<num>"
                                                              "  \tTake one numver - the barrier of contig len, def = 500. Only for ref graph builder" },
            { DIST_BETWEEN_PAIR_READS, 0, "d", "pair_read_dist", Arg::Numeric, " -d <num>, \t --pair_read_dist=<num>"
                                                                                       "\t Take one number - minimal posible len between pair reads"},
            { REFFILE, 0, "r", "ref_file_name", Arg::Required, " -r <arg>, \t --ref_file_name=<arg>"
                                                                       "\t Must have an argument, file name of referense file."},
            { QUERYFILE, 0, "q", "query_file_name", Arg::Required, " -q <arg>, \t --query_file_name=<arg>"
                                                                           "\t Must have an argument, "
                                                                           "file name of contigs will be aligned on ref contigs."
                                                                           "ONLY REF setting"},
            { READSFILE, 0, "p", "read_file",    Arg::Required, " -p <arg>, \t --read_file=<arg> "
                                                                        "\t Must have an argument, file name of reads file."},
            { LIBNAME, 0, "l", "lib_name", Arg::Required, " -l <arg>, \t --lib_name=<arg>"
                                                                  "\t Must have an argument, the name of lib will write on edges"},
            { TSVFILE, 0, "t", "tsv_file_name", Arg::Required, " -t <arg>, \t --tsv_file_name=<arg>"
                                                                       "\t Must have an argument, file name of alignment contigs"
                                                                       "on ref tsv file. ONLY REF setting"},
            { OUTPUT, 0, "o", "wirte_edge_info", Arg::None, "-o, \t --write_edge_info"
                                                                    "\t write add info about edges"},
            { 0, 0, 0, 0, 0, 0 }
    };

public:
    void evaluate(int argc, char **argv); //build conection with this args
};


#endif //SCAFFOLDER_GRAPHCONTROL_H
