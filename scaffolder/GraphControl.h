//
// Created by olga on 08.10.16.
//

#ifndef SCAFFOLDER_GRAPHCONTROL_H
#define SCAFFOLDER_GRAPHCONTROL_H

#include "ConigGraph.h"
#include "GraphBuilder/DNAPairReadGraphBuilder.h"
#include "OutsideLib/optionparser.h"

class GraphControl {
private:
    ConigGraph graph;

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
            string arg = option.arg;
            if (arg == "RNA_PAIR" || arg == "RNA_SPLIT" || arg == "DNA_PAIR") {
                return option::ARG_OK;
            }

            return option::ARG_ILLEGAL;
        }
    };

    enum  optionIndex { UNKNOWN, HELP, NEW, SAM1, SAM2, MINCONTIGLEN, MINEDGELEN};
    const option::Descriptor usage[8] = {
            { UNKNOWN, 0,"", "",        Arg::Unknown, "USAGE: scaffolder [options]\n\n"
                                                              "Options:" },
            { HELP,         0,"", "help",         Arg::None,    "  \t--help  \tPrint usage and exit." },
            { NEW,          0,"n","new",          Arg::NewLib,  "  -n <RNA_PAIR|RNA_SPLIT|DNA_PAIR>, \t--new=<RNA_PAIR|RNA_SPLIT|DNA_PAIR>"
                                                              "  \tNeed to take one of string RNA_PAIR, RNA_SPLIT, DNA_PAIR" },
            { SAM1,         0,"f","sam1",         Arg::Required,"  -f <arg>, \t--sam1=<arg>  \tMust have an argument,"
                                                              "file name of first sam file. def = read1.sam" },
            { SAM2,         0,"s","sam2",         Arg::Required,"  -s <arg>, \t--sam2=<arg>  \tMust have an argument,"
                                                              "file name of second sam file. def = read2.sam" },
            { MINCONTIGLEN, 0,"c","mincontiglen", Arg::Numeric, "  -c <num>, \t--mincontiglen=<num>"
                                                              "  \tTake one numver - the barrier of contig len, def = 0" },
            { MINEDGELEN,   0,"e","minedgelen",   Arg::Numeric, "  -e <num>, \t---minedgelen=<num>"
                                                                        "\t Take one numver - the barrier of edge weight, def = 0"},
            { 0, 0, 0, 0, 0, 0 }
    };


public:
    void evaluate(int argc, char **argv);
};


#endif //SCAFFOLDER_GRAPHCONTROL_H
