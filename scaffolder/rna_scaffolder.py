#!/usr/bin/python3
import sys
import os
import argparse
from shutil import copyfile

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

path_to_exec_dir = os.path.dirname(os.path.abspath(__file__)) + "/"

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--contigs", "-c", nargs=1, dest="contigs", help="path to contigs", type=str)
    parser.add_argument("--rna-p", dest="rnap", nargs=2, help="path to rna pair reads file", type=str, action='append')
    parser.add_argument("--rna-s", dest="rnas", nargs=1, help="path to rna read file", type=str, action='append')
    parser.add_argument("--local_output_dir", "-o", nargs=1, help="use this output dir", type=str)
    args = parser.parse_args()
    return args

def alig_reads(contig_file_name, rnap_list, rnas_list):
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

    for i in range(len(rnap_list)):
        prevdir = os.getcwd();
        lib_name = "rnap" + str(i)
        log.log("START ALIG: " + lib_name)
        lib_dir = os.path.dirname(os.path.abspath(lib_name) + "/")
        if not os.path.exists(lib_dir):
            os.makedirs(lib_dir)
        os.chdir(lib_dir)

        os.system("STAR --runThreadN 20 --genomeDir ../genomeDir --readFilesIn " + rnap_list[i][0])
        os.system("mv Aligned.out.sam rna1.sam")
        os.system("STAR --runThreadN 20 --genomeDir ../genomeDir --readFilesIn " + rnap_list[i][1])
        os.system("mv Aligned.out.sam rna2.sam")
        os.chdir(prevdir)

    for i in range(len(rnas_list)):
        lib_name = "rnas" + str(i)
        lib_dir = os.path.dirname(os.path.abspath(lib_name) + "/")
        if not os.path.exists(lib_dir):
            os.makedirs(lib_dir)

    return

def build_graph(contig_file_name, rnap_list, rnas_list):
    for i in range(len(rnap_list)):
        prevdir = os.getcwd();
        lib_name = "rnap" + str(i)
        log.log("START BUILD GRAPH: " + lib_name)
        lib_dir = os.path.dirname(os.path.abspath(lib_name) + "/")
        os.chdir(lib_dir)

        os.system(path_to_exec_dir + "build -n RNA_SPLIT -r " +
                  contig_file_name + " -p " + rnap_list[i][0] + " -l " + lib_name + "_sp1 " +
                  " -n RNA_SPLIT -r " + contig_file_name + " -p " + rnap_list[i][1] + " -l " + lib_name + "_sp2 " +
                  " -n RNA_PAIR -f rna1.sam -s rna2.sam -l " + lib_name)
        os.chdir(prevdir)

    for i in range(len(rnas_list)):
        prevdir = os.getcwd();
        lib_name = "rnas" + str(i)
        log.log("START BUILD GRAPH: " + lib_name)
        lib_dir = os.path.dirname(os.path.abspath(lib_name) + "/")
        os.chdir(lib_dir)

        os.system(path_to_exec_dir + "build -n RNA_SPLIT -r " + contig_file_name + " -p " + rnap_list[i][0] + " -l " + lib_name)
        os.chdir(prevdir)

    return

def merge_graph(rnap_list, rnas_list):
    args = ""
    for i in range(len(rnap_list)):
        args += "rnap" + str(i) + "/graph.gr "

    for i in range(len(rnas_list)):
        args += "rnas" + str(i) + "/graph.gr "

    args += "graph.gr"
    os.system(path_to_exec_dir + "mergeGraph " + args)

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

    os.system(path_to_exec_dir + "filter " + os.path.abspath("filter_config"))
    return

def create_scaffolds(contig_file_name):
    f = open("filter_config", 'w')
    f.write("uploadGraph out.gr\n")
    f.write("minContig 500\n")
    f.write("mergeSimplePath " + contig_file_name + " scaffolds.fa\n")
    f.write("exit\n")
    f.close()

    os.system(path_to_exec_dir + "filter " + os.path.abspath("filter_config"))
    return

def run(args):
    if args.contigs == None:
        log.err("none contig file provide")
        return

    contig_file_name = os.path.abspath(args.contigs[0])
    rnap_list = []
    rnas_list = []
    if args.rnap:
        for i in range(len(args.rnap)):
            rnap_list.append([os.path.abspath(args.rnap[i][0]), os.path.abspath(args.rnap[i][1])])

    if args.rnas:
        for i in range(len(args.rnas)):
            rnas_list.append(os.path.abspath(args.rnas[i][0]))


    main_out_dir = os.path.abspath(args.local_output_dir[0]) + "/"

    if args.local_output_dir != None:
        out_dir = main_out_dir + "tmp/"
        log.log("OUTPUT DIR: " + out_dir)
        directory = os.path.dirname(out_dir)
        if not os.path.exists(directory):
            log.log("MKDIR")
            os.makedirs(directory)
        os.chdir(directory)

    alig_reads(contig_file_name, rnap_list, rnas_list)
    build_graph(contig_file_name, rnap_list, rnas_list)
    merge_graph(rnap_list, rnas_list)
    create_scaffolds(contig_file_name)

    directory = os.path.dirname(main_out_dir)
    os.chdir(directory)

    copyfile("tmp/scaffolds.fa", "scaffolds.fa")
    copyfile("tmp/out.gr", "graph.gr")
    copyfile("tmp/out.info", "out.info")
    return

args = parse_args()
run(args)