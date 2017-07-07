#!/usr/bin/python3
import sys
import os
import argparse

#Log class, use it, not print
class Log:
    text = ""

    def log(self, s):
        self.text += s + "\n"
        print(s)

    def warn(self, s):
        msg = "WARNING: " + s
        self.text += msg + "\n"
        sys.stdout.write(msg)
        sys.stdout.flush()

    def err(self, s):
        msg = "ERROR: " + s + "\n"
        self.text += msg
        sys.stdout.write(msg)
        sys.stdout.flush()

    def print_log(self):
        print(self.text)

    def get_log(self):
        return self.text

log = Log()

class Lib:
    def __init__(self, path, id, name):
        self.path = []
        for p in path:
            self.path.append(os.path.abspath(p))
        self.id = id
        self.color = "#000000"
        self.label = name + "_" + str(id)
        self.name = name + "_" + str(id)

libsType = {"rnap", "rnas", "dnap", "ref", "selfinfo", "selfpath"}

class StoreArgAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if not 'lib_cnt' in namespace:
            setattr(namespace, 'lib_cnt', 0)
        lib_cnt = namespace.lib_cnt

        if not 'libs' in namespace:
            libs = dict()
            for lib in libsType:
                libs[lib] = []

            setattr(namespace, 'libs', libs)

        libs = namespace.libs
        libs[self.dest].append(Lib(values, lib_cnt, self.dest))

        lib_cnt += 1
        setattr(namespace, 'libs', libs)
        setattr(namespace, 'lib_cnt', lib_cnt)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--contigs", "-c", nargs=1, dest="contigs", help="path to contigs", type=str)
    parser.add_argument("--rna-p", dest="rnap", nargs=2, help="path to rna pair reads file", type=str, action=StoreArgAction)
    parser.add_argument("--rna-s", dest="rnas", nargs=1, help="path to rna read file", type=str, action=StoreArgAction)
    parser.add_argument("--local_output_dir", "-o", nargs=1, help="use this output dir", type=str)
    parser.add_argument("--ref", dest="ref", nargs=1, help="path to reference", type=str, action=StoreArgAction)
    parser.add_argument("--dna-p", dest="dnap", nargs=2, help="path to dna pair reads file", type=str, action=StoreArgAction)
    parser.add_argument("--scafinfo", nargs=1, help="path to .info file with info about scaffolds", type=str, action=StoreArgAction)
    parser.add_argument("--scafpath", nargs=1, help="path to both.path file with info about scaffolds", type=str, action=StoreArgAction)
    parser.add_argument("--label", "-l", nargs='*', help="list with labels for all libs in definition order", type=str)
    parser.add_argument("--color", nargs='*', help="list with color for all libs in definition order", type=str)
    args = parser.parse_args()
    return args

def alig_reads(contig_file_name, args):
    genome_dir = "genomeDir"
    gen_dir = os.path.dirname(os.path.abspath(genome_dir) + "/")
    if not os.path.exists(gen_dir):
        os.makedirs(gen_dir)

    try:
        os.system("STAR --runMode genomeGenerate --genomeDir genomeDir --runThreadN 20 --genomeSAindexNbases 10 --genomeFastaFiles " +
                  contig_file_name + " --limitGenomeGenerateRAM 90000000000")
    except:
        log.err(sys.exc_info()[0])
        return

    for i in range(len(args.libs["rnap"])):
        prevdir = os.getcwd();
        log.log("START ALIG: " + args.libs["rnap"][i].label)
        lib_dir = os.path.dirname(os.path.abspath(args.libs["rnap"][i].name) + "/")
        os.chdir(lib_dir)

        os.system("STAR --runThreadN 20 --genomeDir ../genomeDir --readFilesIn " + args.libs["rnap"][i].path[0])
        os.system("mv Aligned.out.sam rna1.sam")
        os.system("STAR --runThreadN 20 --genomeDir ../genomeDir --readFilesIn " + args.libs["rnap"][i].path[1])
        os.system("mv Aligned.out.sam rna2.sam")
        os.chdir(prevdir)

    for i in range(len(args.libs["dnap"])):
        prevdir = os.getcwd();
        log.log("START ALIG: " + args.libs["dnap"][i].label)
        lib_dir = os.path.dirname(os.path.abspath(args.libs["dnap"][i].name) + "/")
        os.chdir(lib_dir)

        os.system("bowtie2-build " + contig_file_name + " contig")
        if args.libs["dnap"][i].path[0].endswith(".fa") or args.libs["dnap"][i].path[0].endswith(".fasta") or args.libs["dnap"][i].path[0].endswith(".mfa") or args.libs["dnap"][i].path[0].endswith(".fna"):
            os.system("bowtie2 -x contig -f --ignore-quals -U " + args.libs["dnap"][i].path[0] + " -S dna1.sam")
            os.system("bowtie2 -x contig -f --ignore-quals -U " + args.libs["dnap"][i].path[1] + " -S dna2.sam")
        else:
            os.system("bowtie2 -x contig -U " + args.libs["dnap"][i].path[0] + " -S dna1.sam")
            os.system("bowtie2 -x contig -U " + args.libs["dnap"][i].path[1] + " -S dna2.sam")
        os.chdir(prevdir)

    return

def build_graph(contig_file_name, args):
    for lib in args.libs["rnap"]:
        prevdir = os.getcwd();
        log.log("START BUILD GRAPH: " + lib.label)
        lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
        os.chdir(lib_dir)

        os.system("../../build -n RNA_SPLIT -r " +
                  contig_file_name + " -p " + lib.path[0] + " -l " + lib.label + "_sp1 " +
                  " -n RNA_SPLIT -r " + contig_file_name + " -p " + lib.path[1] + " -l " + lib.label + "_sp2 " +
                  " -n RNA_PAIR -f rna1.sam -s rna2.sam -l " + lib.label)
        os.chdir(prevdir)

    for lib in args.libs["rnas"]:
        prevdir = os.getcwd();
        log.log("START BUILD GRAPH: " + lib.label)
        lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
        os.chdir(lib_dir)
        os.system("../../build -n RNA_SPLIT -r " + contig_file_name + " -p " + lib.path[0] + " -l " + lib.label)
        os.chdir(prevdir)

    for lib in args.libs["dnap"]:
        prevdir = os.getcwd();
        log.log("START BUILD GRAPH: " + lib.label)
        lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
        os.chdir(lib_dir)
        os.system("../../build -n DNA_PAIR -f dna1.sam -s dna2.sam -l " + lib.label)
        os.chdir(prevdir)

    for lib in args.libs["ref"]:
        prevdir = os.getcwd();
        log.log("START BUILD GRAPH: " + lib.label)
        lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
        os.chdir(lib_dir)
        os.system("../../build -n REF -r " + lib.path[0] + " -q " + contig_file_name + " -l " + lib.label)
        os.chdir(prevdir)

    return

def merge_graph(rnap_list, rnas_list):
    args = ""
    for i in range(len(rnap_list)):
        args += "rnap" + str(i) + "/graph.gr "

    for i in range(len(rnas_list)):
        args += "rnas" + str(i) + "/graph.gr "

    args += "graph.gr"
    os.system("../mergeGraph " + args)

    list_s50 = []
    list_s30 = []
    list_p = []

    cur_lib = 0
    for i in range(len(rnap_list)):
        list_s50.append(cur_lib)
        cur_lib += 1
        list_s30.append(cur_lib)
        cur_lib += 1
        list_s50.append(cur_lib)
        cur_lib += 1
        list_s30.append(cur_lib)
        cur_lib += 1
        list_p.append(cur_lib)
        cur_lib += 1

    for i in range(len(rnas_list)):
        list_s50.append(cur_lib)
        cur_lib += 1
        list_s30.append(cur_lib)
        cur_lib += 1


    f = open("filter_config", 'w')
    f.write("uploadGraph graph.gr\n")
    for i in range(len(list_s50)-1, 0, -1):
        f.write("mergeLib " + str(list_s50[i]) + " " + str(list_s50[i - 1]) + " sp_50\n")
    for i in range(len(list_s30)-1, 0, -1):
        f.write("mergeLib " + str(list_s30[i]) + " " + str(list_s30[i - 1]) + " sp_30\n")
    for i in range(len(list_p)-1, 0, -1):
        f.write("mergeLib " + str(list_p[i]) + " " + str(list_p[i - 1]) + " pair\n")

    f.write("print out.gr\n")
    f.write("exit\n")
    f.close()

    os.system("../filter " + os.path.abspath("filter_config"))
    return

def run(args):
    if args.contigs == None:
        log.err("none contig file provide")
        return

    contig_file_name = os.path.abspath(args.contigs[0])

    if args.local_output_dir != None:
        out_dir = os.path.abspath(args.local_output_dir[0]) + "/"
        log.log("OUTPUT DIR: " + out_dir)
        directory = os.path.dirname(out_dir)
        if not os.path.exists(directory):
            log.log("MKDIR")
            os.makedirs(directory)
        os.chdir(directory)

    print(args.libs)

    if args.color != None and len(args.color) != args.lib_cnt:
        log.err("wrong number of color provide")

    if args.label != None and len(args.label) != args.lib_cnt:
        log.err("wrong number of labels provide")

    for lib_type in libsType:
        for lib in args.libs[lib_type]:
            lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
            if not os.path.exists(lib_dir):
                os.makedirs(lib_dir)

            if args.color != None:
                lib.color = args.color[lib.id]
            if args.label != None:
                lib.label = args.label[lib.id]


    alig_reads(contig_file_name, args)
    build_graph(contig_file_name, args)
    #merge_graph(args)


    return

args = parse_args()
run(args)