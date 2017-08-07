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

def alig_split(lib_name, reads, flag):
    prevdir = os.getcwd()
    log.log("START ALIG: " + lib_name)
    lib_dir = os.path.dirname(os.path.abspath(lib_name) + "/")
    if not os.path.exists(lib_dir):
        os.makedirs(lib_dir)
    os.chdir(lib_dir)

    os.system(path_to_exec_dir + "readSplitter " + str(flag) + " " + reads + " reads1.fasta reads2.fasta")
    os.system("STAR --runThreadN 20 --genomeDir ../genomeDir --readFilesIn reads1.fasta")
    os.system("mv Aligned.out.sam rna1.sam")

    os.system("STAR --runThreadN 20 --genomeDir ../genomeDir --readFilesIn reads2.fasta")
    os.system("mv Aligned.out.sam rna2.sam")

    os.chdir(prevdir)

def alig_pair_reads(i, rnap1, rnap2):
    prevdir = os.getcwd()
    lib_name = "rnap" + str(i)
    log.log("START ALIG: " + lib_name)
    lib_dir = os.path.dirname(os.path.abspath(lib_name) + "/")
    if not os.path.exists(lib_dir):
        os.makedirs(lib_dir)
    os.chdir(lib_dir)

    unm1 = ""
    unm2 = ""

    os.system("STAR --runThreadN 20 --genomeDir ../genomeDir --outReadsUnmapped Fastx --readFilesIn " + rnap1)
    os.system("mv Aligned.out.sam rna1.sam")
    if rnap1[-1] == "q":
        os.system("mv Unmapped.out.mate1 Unmapped1.fastq")
        unm1 = "../" + lib_name + "/Unmapped1.fastq"
    else:
        os.system("mv Unmapped.out.mate1 Unmapped1.fasta")
        unm1 = "../" + lib_name + "/Unmapped1.fasta"

    os.system("STAR --runThreadN 20 --genomeDir ../genomeDir --outReadsUnmapped Fastx --readFilesIn " + rnap2)
    os.system("mv Aligned.out.sam rna2.sam")
    if rnap1[-1] == "q":
        os.system("mv Unmapped.out.mate1 Unmapped2.fastq")
        unm2 = "../" + lib_name + "/Unmapped2.fastq"
    else:
        os.system("mv Unmapped.out.mate1 Unmapped2.fasta")
        unm2 = "../" + lib_name + "/Unmapped2.fasta"

    os.chdir(prevdir)

    alig_split("rnap" + str(i) + "_50_1", unm1, 0)
    alig_split("rnap" + str(i) + "_50_2", unm2, 0)
    alig_split("rnap" + str(i) + "_30_1", "../" + lib_name + "/rna1.sam", 1)
    alig_split("rnap" + str(i) + "_30_2", "../" + lib_name + "/rna2.sam", 1)


def alig_single_reads(i, rnas):
    prevdir = os.getcwd()
    lib_name = "rnas" + str(i) + "_50"
    lib_dir = os.path.dirname(os.path.abspath(lib_name) + "/")
    if not os.path.exists(lib_dir):
        os.makedirs(lib_dir)
    os.chdir(lib_name)
    unm = ""
    os.system("STAR --runThreadN 20 --genomeDir ../genomeDir --outReadsUnmapped Fastx --readFilesIn " + rnas)
    os.system("mv Aligned.out.sam rna.sam")
    if rnas[-1] == "q":
        os.system("mv Unmapped.out.mate1 Unmapped.fastq")
        unm = "Unmapped.fastq"
    else:
        os.system("mv Unmapped.out.mate1 Unmapped.fasta")
        unm = "Unmapped.fasta"


    os.chdir(prevdir)
    alig_split(lib_name, unm, 0)
    alig_split("rnas" + str(i) + "_30", "../rnas" + str(i) + "_30/rna.sam", 1)

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
        alig_pair_reads(i, rnap_list[i][0], rnap_list[i][1])


    for i in range(len(rnas_list)):
        alig_single_reads(i, rnas_list[i][0])
    return

def runGraphBuilder(lib_name, prevdir, type):
    log.log("START BUILD GRAPH: " + lib_name)
    lib_dir = os.path.dirname(os.path.abspath(lib_name) + "/")
    os.chdir(lib_dir)
    os.system(path_to_exec_dir + "build " + type + " rna1.sam rna2.sam " + lib_name)
    os.chdir(prevdir)
    return


def build_graph(contig_file_name, rnap_list, rnas_list):
    for i in range(len(rnap_list)):
        prevdir = os.getcwd();
        runGraphBuilder("rnap" + str(i), prevdir, "RNA_PAIR")
        runGraphBuilder("rnap" + str(i) + "_50_1", prevdir, "RNA_SPLIT_50")
        runGraphBuilder("rnap" + str(i) + "_50_2", prevdir, "RNA_SPLIT_50")
        runGraphBuilder("rnap" + str(i) + "_30_1", prevdir, "RNA_SPLIT_30")
        runGraphBuilder("rnap" + str(i) + "_30_2", prevdir, "RNA_SPLIT_30")

    for i in range(len(rnas_list)):
        prevdir = os.getcwd();
        runGraphBuilder("rnas" + str(i) + "_50", prevdir, "RNA_SPLIT_50")
        runGraphBuilder("rnas" + str(i) + "_30", prevdir, "RNA_SPLIT_30")
    return

def merge_graph(rnap_list, rnas_list):
    args = ""
    for i in range(len(rnap_list)):
        args += "rnap" + str(i) + "_50_1/graph.gr "
        args += "rnap" + str(i) + "_30_1/graph.gr "
        args += "rnap" + str(i) + "_50_2/graph.gr "
        args += "rnap" + str(i) + "_30_2/graph.gr "
        args += "rnap" + str(i) + "/graph.gr "


    for i in range(len(rnas_list)):
        args += "rnas" + str(i) + "_50/graph.gr "
        args += "rnas" + str(i) + "_30/graph.gr "

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


    main_out_dir = os.path.abspath(".") + "/"
    if args.local_output_dir != None:
        main_out_dir = os.path.abspath(args.local_output_dir[0]) + "/"

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