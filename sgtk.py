#!/usr/bin/env python3
import sys
import os
import argparse
import json
from shutil import copyfile
from shutil import move
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

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
        parser.print_help()
        sys.exit()

    def print_log(self):
        print(self.text)

    def get_log(self):
        return self.text

log = Log()
path_to_exec_dir = os.path.dirname(os.path.abspath(__file__)) + "/"

class Lib:
    def __init__(self, path, id, name):
        self.path = []
        for p in path:
            self.path.append(os.path.abspath(p))
        self.id = id
        self.color = "#000000"
        self.label = name + "_" + str(id)
        self.name = name + "_" + str(id)

    def fix_graph_file(self):
        if self.name.startswith("scaf"):
            return

        copyfile(self.name + "/graph.gr", "tmp")

        g = open(self.name + "/graph.gr", "w")
        f = open("tmp", "r")

        if not self.name.startswith("rna"):
            f.readline()
            g.write("1\n")
            libinfo = f.readline().split(' ')
            libinfo[2] = self.color
            libinfo[3] = self.label
            g.write(' '.join(libinfo))
            str = f.read()
            g.write(str)
        elif self.name.startswith("rnap"):
            f.readline()
            g.write("3\n")
            libinfo = f.readline().split(' ')
            libinfo[2] = self.color
            libinfo[3] = self.label + "_sp50"
            g.write(' '.join(libinfo))

            libinfo = f.readline().split(' ')
            libinfo[2] = self.color
            libinfo[3] = self.label + "_sp30"
            g.write(' '.join(libinfo))

            libinfo = f.readline().split(' ')
            libinfo[2] = self.color
            libinfo[3] = self.label + "_pair"
            g.write(' '.join(libinfo))

            str = f.read()
            g.write(str)
        else:
            f.readline()
            g.write("2\n")
            libinfo = f.readline().split(' ')
            libinfo[2] = self.color
            libinfo[3] = self.label + "_sp50"
            g.write(' '.join(libinfo))

            libinfo = f.readline().split(' ')
            libinfo[2] = self.color
            libinfo[3] = self.label + "_sp30"
            g.write(' '.join(libinfo))

            str = f.read()
            g.write(str)

        f.close()
        g.close()

libsType = {"rnap", "rnas", "rf", "ff", "scg", "ref", "scafinfo", "scafpath", "scaffolds", "refcoord", "fr", "long", "fastg", "gfa", "frsam", "rfsam", "ffsam"}

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
    global parser
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--local_output_dir", "-o", nargs=1, help="output folder", type=str)

    parser.add_argument("--contigs", "-c", nargs=1, dest="contigs", help="path to contigs  in FASTA format", type=str, action='append')
    parser.add_argument("--scaffolds", "-s", nargs=1, dest="scaffolds", help="path to scaffolds in FASTA format", type=str, action=StoreArgAction)
    parser.add_argument("--fastg", nargs=1, dest="fastg", help="path to assembly graph in FASTG format", type=str, action=StoreArgAction)
    parser.add_argument("--gfa", nargs=1, dest="gfa", help="path to assembly graph in GFA1 format", type=str, action=StoreArgAction)

    parser.add_argument("--fr", dest="fr", nargs=2, help="paths to forward-reverse paired reads in FASTQ/FASTA format", type=str, action=StoreArgAction)
    parser.add_argument("--rf", dest="rf", nargs=2, help="paths to reverse-forward paired reads in FASTQ/FASTA format", type=str, action=StoreArgAction)
    parser.add_argument("--ff", dest="ff", nargs=2, help="paths to forward-forward paired reads in FASTQ/FASTA format", type=str, action=StoreArgAction)

    parser.add_argument("--fr_sam", dest="frsam", nargs=2, help="paths to alignment of forward-reverse paired reads in SAM/BAM format", type=str, action=StoreArgAction)
    parser.add_argument("--rf_sam", dest="rfsam", nargs=2, help="paths to alignment of reverse-forward paired reads in SAM/BAM format", type=str, action=StoreArgAction)
    parser.add_argument("--ff_sam", dest="ffsam", nargs=2, help="paths to alignment of forward-forward paired reads in SAM/BAM format", type=str, action=StoreArgAction)


    parser.add_argument("--long", dest="long", nargs=1, help="path to long reads (PacBio/Oxford Nanopore) file in FASTQ/FASTA format", type=str, action=StoreArgAction)
    parser.add_argument("--rna-p", dest="rnap", nargs=2, help="paths to RNA-Seq paired reads  in SAM/BAM format", type=str, action=StoreArgAction)
    parser.add_argument("--rna-s", dest="rnas", nargs=1, help="path to RNA-Seq single reads  in SAM/BAM format", type=str, action=StoreArgAction)

    parser.add_argument("--ref", dest="ref", nargs=1, help="path to reference genome in FASTA format", type=str, action='append')
    parser.add_argument("--refcoord", dest="refcoord", nargs=2, help="path to ref and to alignment of contigs to reference in coord format", type=str, action='append')

    parser.add_argument("--gr", nargs=1, dest="graph", help="path to graph in .gr format", type=str, action='append')
    parser.add_argument("--scg", nargs=1, dest="scg", help="path to file with connection list", type=str, action=StoreArgAction)
    parser.add_argument("--scafinfo", nargs=1, help="path to .info file with information about scaffolds", type=str, action=StoreArgAction)
    parser.add_argument("--scafpath", nargs=1, help="path to .paths file with information about scaffolds", type=str, action=StoreArgAction)

    parser.add_argument("--label", "-l", nargs='*', help="list with labels for all sorces in the corresponding order", type=str, action='store')
    parser.add_argument("--color", nargs='*', help="list with colors for all  sorces in the corresponding order", type=str, action='store')
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

def alig_pair_rna_reads(rnap):
    prevdir = os.getcwd()
    log.log("START ALIG: " + rnap.label)
    lib_dir = os.path.dirname(os.path.abspath(rnap.name) + "/")
    os.chdir(lib_dir)

    unm1 = ""
    unm2 = ""

    os.system("STAR --runThreadN 4 --genomeDir ../genomeDir --outReadsUnmapped Fastx --readFilesIn " + rnap.path[0])
    os.system("mv Aligned.out.sam rna1.sam")
    if rnap.path[0][-1] == "q":
        os.system("mv Unmapped.out.mate1 Unmapped1.fastq")
        unm1 = "../" + rnap.name + "/Unmapped1.fastq"
    else:
        os.system("mv Unmapped.out.mate1 Unmapped1.fasta")
        unm1 = "../" + rnap.name + "/Unmapped1.fasta"

    os.system("STAR --runThreadN 4 --genomeDir ../genomeDir --outReadsUnmapped Fastx --readFilesIn " + rnap.path[1])
    os.system("mv Aligned.out.sam rna2.sam")

    if rnap.path[1][-1] == "q":
        os.system("mv Unmapped.out.mate1 Unmapped2.fastq")
        unm2 = "../" + rnap.name + "/Unmapped2.fastq"
    else:
        os.system("mv Unmapped.out.mate1 Unmapped2.fasta")
        unm2 = "../" + rnap.name + "/Unmapped2.fasta"

    os.chdir(prevdir)

    alig_split(rnap.name + "_50_1", unm1, 0)
    alig_split(rnap.name + "_50_2", unm2, 0)
    alig_split(rnap.name + "_30_1", "../" + rnap.name + "/rna1.sam", 1)
    alig_split(rnap.name + "_30_2", "../" + rnap.name + "/rna2.sam", 1)

def alig_pair_dna_reads(dnap, contig_file_name):
    prevdir = os.getcwd()
    log.log("START ALIG: " + dnap.label)
    lib_dir = os.path.dirname(os.path.abspath(dnap.name) + "/")
    os.chdir(lib_dir)

    os.system("minimap2 -t 16 -ax sr " + contig_file_name + " " + dnap.path[0] + " > dna1.sam")
    os.system("minimap2 -t 16 -ax sr " + contig_file_name + " " + dnap.path[1] + " > dna2.sam")

    #os.system("bowtie2-build " + contig_file_name + " contig")
    #if dnap.path[0].endswith(".fa") or dnap.path[0].endswith(".fasta") or dnap.path[0].endswith(".mfa") or dnap.path[0].endswith(".fna"):
    #    os.system("bowtie2 -x contig -f --ignore-quals -U " + dnap.path[0] + " -S dna1.sam")
    #    os.system("bowtie2 -x contig -f --ignore-quals -U " + dnap.path[1] + " -S dna2.sam")
    #else:
    #    os.system("bowtie2 -x contig -U " + dnap.path[0] + " -S dna1.sam")
    #    os.system("bowtie2 -x contig -U " + dnap.path[1] + " -S dna2.sam")
    os.chdir(prevdir)

def alig_long(long, contig_file_name):
    prevdir = os.getcwd()
    log.log("START ALIG: " + long.label)
    lib_dir = os.path.dirname(os.path.abspath(long.name) + "/")
    os.chdir(lib_dir)

    os.system("minimap2 -t 16 -x map-pb " + contig_file_name + " " + long.path[0] + " > out.paf")

    os.chdir(prevdir)

def alig_single_rna_reads(rnas):
    prevdir = os.getcwd()
    lib_name = rnas.name + "_50"
    lib_dir = os.path.dirname(os.path.abspath(lib_name) + "/")
    if not os.path.exists(lib_dir):
        os.makedirs(lib_dir)
    os.chdir(lib_name)
    unm = ""
    os.system("STAR --runThreadN 4 --genomeDir ../genomeDir --outReadsUnmapped Fastx --readFilesIn " + rnas)
    os.system("mv Aligned.out.sam rna.sam")
    if rnas.path[0][-1] == "q":
        os.system("mv Unmapped.out.mate1 Unmapped.fastq")
        unm = "Unmapped.fastq"
    else:
        os.system("mv Unmapped.out.mate1 Unmapped.fasta")
        unm = "Unmapped.fasta"

    os.chdir(prevdir)
    alig_split(lib_name, unm, 0)
    alig_split(rnas.name + "_30", "../" + rnas.name + "_30/rna.sam", 1)
    return

def alig_reads(contig_file_name, args):
    genome_dir = "genomeDir"
    gen_dir = os.path.dirname(os.path.abspath(genome_dir) + "/")
    if not os.path.exists(gen_dir):
        os.makedirs(gen_dir)

    try:
        if (len(args.libs["rnap"]) > 0 or len(args.libs["rnas"]) > 0):
            log.log("STAR --runMode genomeGenerate --genomeDir genomeDir --runThreadN 20 --genomeSAindexNbases 10 --genomeFastaFiles " +
                    contig_file_name + " --limitGenomeGenerateRAM 90000000000 --genomeChrBinNbits 15")
            os.system("STAR --runMode genomeGenerate --genomeDir genomeDir --runThreadN 20 --genomeSAindexNbases 10 --genomeFastaFiles " +
                      contig_file_name + " --limitGenomeGenerateRAM 90000000000 --genomeChrBinNbits 15")
    except:
        log.err(sys.exc_info()[0])
        return

    for i in range(len(args.libs["rnap"])):
        alig_pair_rna_reads(args.libs["rnap"][i])

    for i in range(len(args.libs["fr"])):
        alig_pair_dna_reads(args.libs["fr"][i], contig_file_name)

    for i in range(len(args.libs["rf"])):
        alig_pair_dna_reads(args.libs["rf"][i], contig_file_name)

    for i in range(len(args.libs["ff"])):
        alig_pair_dna_reads(args.libs["ff"][i], contig_file_name)

    for i in range(len(args.libs["long"])):
        alig_long(args.libs["long"][i], contig_file_name)

    for i in range(len(args.libs["rnas"])):
        alig_single_rna_reads(args.libs["rnas"][i])

    return

def runGraphBuilder(lib_name, prevdir, type, label):
    log.log("START BUILD GRAPH: " + lib_name)
    lib_dir = os.path.dirname(os.path.abspath(lib_name) + "/")
    os.chdir(lib_dir)
    os.system(path_to_exec_dir + "buildApp " + type + " rna1.sam rna2.sam " + label)
    os.chdir(prevdir)
    return

def build_graph(contig_file_name, args):
    for lib in args.libs["rnap"]:
        prevdir = os.getcwd()
        runGraphBuilder(lib.name, prevdir, "RNA_PAIR", lib.label)
        runGraphBuilder(lib.name + "_50_1", prevdir, "RNA_SPLIT_50", lib.label)
        runGraphBuilder(lib.name + "_50_2", prevdir, "RNA_SPLIT_50", lib.label)
        runGraphBuilder(lib.name + "_30_1", prevdir, "RNA_SPLIT_30", lib.label)
        runGraphBuilder(lib.name + "_30_2", prevdir, "RNA_SPLIT_30", lib.label)

    for lib in args.libs["rnas"]:
        prevdir = os.getcwd()
        runGraphBuilder(lib.name + "_50", prevdir, "RNA_SPLIT_50", lib.label)
        runGraphBuilder(lib.name + "_30", prevdir, "RNA_SPLIT_30", lib.label)

    for lib in args.libs["fr"]:
        prevdir = os.getcwd()
        log.log("START BUILD GRAPH: " + lib.label)
        lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
        os.chdir(lib_dir)
        os.system(path_to_exec_dir + "buildApp DNA_PAIR_FR dna1.sam dna2.sam " + lib.label)
        os.chdir(prevdir)

    for lib in args.libs["rf"]:
        prevdir = os.getcwd()
        log.log("START BUILD GRAPH: " + lib.label)
        lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
        os.chdir(lib_dir)
        os.system(path_to_exec_dir + "buildApp DNA_PAIR_RF dna1.sam dna2.sam " + lib.label)
        os.chdir(prevdir)

    for lib in args.libs["ff"]:
        prevdir = os.getcwd()
        log.log("START BUILD GRAPH: " + lib.label)
        lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
        os.chdir(lib_dir)
        os.system(path_to_exec_dir + "buildApp DNA_PAIR_FF dna1.sam dna2.sam " + lib.label)
        os.chdir(prevdir)


    for lib in args.libs["frsam"]:
        prevdir = os.getcwd()
        log.log("START BUILD GRAPH: " + lib.label)
        lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
        os.chdir(lib_dir)
        os.system(path_to_exec_dir + "buildApp DNA_PAIR_FR " + lib.path[0] + " " + lib.path[1]  + " " + lib.label)
        os.chdir(prevdir)

    for lib in args.libs["rfsam"]:
        prevdir = os.getcwd()
        log.log("START BUILD GRAPH: " + lib.label)
        lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
        os.chdir(lib_dir)
        os.system(path_to_exec_dir + "buildApp DNA_PAIR_RF " + lib.path[0] + " " + lib.path[1]  + " " + lib.label)
        os.chdir(prevdir)

    for lib in args.libs["ffsam"]:
        prevdir = os.getcwd()
        log.log("START BUILD GRAPH: " + lib.label)
        lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
        os.chdir(lib_dir)
        os.system(path_to_exec_dir + "buildApp DNA_PAIR_FF " + lib.path[0] + " " + lib.path[1]  + " " + lib.label)
        os.chdir(prevdir)


    for lib in args.libs["long"]:
        prevdir = os.getcwd()
        log.log("START BUILD GRAPH: " + lib.label)
        lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
        os.chdir(lib_dir)
        os.system(path_to_exec_dir + "buildApp LONG out.paf " + contig_file_name + " " + lib.label)
        os.chdir(prevdir)

    for lib in args.libs["scg"]:
        prevdir = os.getcwd()
        log.log("START BUILD GRAPH: " + lib.label)
        lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
        os.chdir(lib_dir)
        os.system(path_to_exec_dir + "buildApp CONNECTION " + lib.path[0] + " " + contig_file_name + " " + lib.label)
        os.chdir(prevdir)

    for lib in args.libs["fastg"]:
        prevdir = os.getcwd()
        log.log("START BUILD GRAPH: " + lib.label)
        lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
        os.chdir(lib_dir)
        os.system(path_to_exec_dir + "buildApp FASTG " + lib.path[0] + " " + contig_file_name + " " + lib.label)
        os.chdir(prevdir)

    for lib in args.libs["gfa"]:
        prevdir = os.getcwd()
        log.log("START BUILD GRAPH: " + lib.label)
        lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
        os.chdir(lib_dir)
        os.system(path_to_exec_dir + "buildApp GFA " + lib.path[0] + " " + lib.label)
        os.chdir(prevdir)
        
    return

def merge_lib_rna(libs):
    f = open("filter_config", 'w')
    f.write("uploadGraph graph.gr\n")
    f.write("mergeLib 2 0 sp_50\n")
    f.write("mergeLib 3 1 sp_30\n")
    f.write("print graph.gr\n")
    f.write("exit\n")
    f.close()

    for lib in libs:
        os.system(path_to_exec_dir + "mergeGraph " +
                  lib.name + "_50_1/graph.gr " +
                  lib.name + "_30_1/graph.gr " +
                  lib.name + "_50_2/graph.gr " +
                  lib.name + "_30_2/graph.gr " +
                  lib.name + "/graph.gr " +
                  lib.name + "/gr.gr")

        copyfile(lib.name + "/gr.gr", lib.name + "/graph.gr")
        prevdir = os.getcwd()
        lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
        os.chdir(lib_dir)
        os.system(path_to_exec_dir + "filterApp " + os.path.abspath("../filter_config"))
        os.chdir(prevdir)
    return

def merge_graph(args):
    if 'libs' in args:
        merge_lib_rna(args.libs["rnap"])

    merge_list = ""

    if 'libs' in args:
        for lib_type in libsType:
            for lib in args.libs[lib_type]:
                if lib_type == "rnap" or lib_type == "rnas" or lib_type == "fr" or lib_type == "rf" or \
                        lib_type == "long" or lib_type == "ff" or lib_type == "scg" or lib_type=="fastg" \
                        or lib_type=="gfa" or lib_type == "frsam" or lib_type == "rfsam" or lib_type == "ffsam":
                    lib.fix_graph_file()
                    if lib_type != "rnas":
                        merge_list += lib.name + "/graph.gr "
                    else:
                        merge_list += lib.name + "_50/graph.gr "
                        merge_list += lib.name + "_30/graph.gr "

    if args.graph != None:
        for gr in args.graph:
            merge_list += gr + " "

    if merge_list == "":
        return
    
    merge_list += "graph.gr"
    os.system(path_to_exec_dir + "mergeGraph " + merge_list)
    return

idbyname = dict()
partnametoname = dict()
lenbyid = []
cntedge = 0
cntlib = 0

output_json = dict()


#save id to idbyname
#write node info to data
def gen_id_from_contig_file(contig_file_name):
    fasta_seq = SeqIO.parse(open(contig_file_name), 'fasta')
    id = 0
    output_json['nodes'] = []
    for fasta in fasta_seq:
        name, lenn = fasta.id, len(fasta.seq.tostring())
        idbyname[name] = id
        idbyname[name + "-rev"] = id + 1
        lenbyid.append(lenn)
        lenbyid.append(lenn)
        output_json['nodes'].append({'id': id, 'name': name, 'len': lenn})
        output_json['nodes'].append({'id': id + 1, 'name': name + "-rev", 'len': lenn})
        id += 2


def save_scaffolds_from_info(lib):
    global cntedge
    global cntlib
    global idbyname
    with open(lib.path[0]) as g:
        output_json['libs'].append({'id': cntlib, 'color': lib.color, 'name': lib.label, 'type': 'SCAFF', 'scaffolds': []})
        scafnum = 0

        for line in g:
            tokens = line.split(" ")
            if (tokens[len(tokens) - 1] == '\n'):
                tokens.pop()

            output_json['libs'][-1]['scaffolds'].append({'name': tokens[0][1:], 'edges': []})
            output_json['libs'][-1]['scaffolds'].append({'name': tokens[0][1:] + "-rev", 'edges': []})

            nodeslist = []
            for i in range(1, len(tokens), 3):
                nm = tokens[i][1:]
                if (tokens[i + 2][0] == '+'):
                    nodeslist.append(idbyname[nm])
                else:
                    nodeslist.append(idbyname[nm]^1)


            for i in range(1, len(nodeslist)):
                output_json['libs'][-1]['scaffolds'][-2]['edges'].append({'id': cntedge, 'from': nodeslist[i - 1],
                                                                          'to': nodeslist[i], 'weight': 1})
                cntedge += 1

            for i in range(len(nodeslist)-2, -1, -1):
                output_json['libs'][-1]['scaffolds'][-1]['edges'].append({'id': cntedge, 'from': nodeslist[i + 1]^1,
                                                                          'to': nodeslist[i]^1, 'weight': 1})
                cntedge += 1

            scafnum += 2
    cntlib += 1


def sortcmp(x, y):
    if (x[0] < y[0]):
        return -1
    if (x[0] == y[0] and x[1] > y[1]):
        return -1
    if (x[0] == y[0] and x[1] == y[1]):
        return 0
    return 1

def save_lens_from_sam(scafflen_by_name, file_name):
    with open(file_name) as f:
        for line in f:
            if (line[0] != '@'):
                return
            if (line[0:3] != '@SQ'):
                continue
            parts = line.split()
            scafname = ""
            curlen = 0
            for part in parts:
                if (part[0:2] == 'LN'):
                    curlen = int(part[3:])
                if (part[0:2] == 'SN'):
                    scafname = part[3:]
            scafflen_by_name[scafname] = curlen
    return

def parse_cigar(cigar):
    letCnt = {'S': 0, 'H': 0, 'M': 0, '=': 0, 'X': 0, 'I': 0, 'D': 0, 'N': 0}
    num = 0
    for i in range(len(cigar)):
        if (cigar[i].isdigit()):
            num = int(num)*10 + int(cigar[i])
        else:
            letCnt[cigar[i]] += num
            num = 0
    return letCnt


def get_align_from_sam_line(line):
    tokens = line.split("\t")
    if (len(tokens) < 10):
        return -1, -1, -1, -1, "", ""

    if (tokens[len(tokens) - 1] == '\n'):
        tokens.pop()
    if (tokens[len(tokens) - 1][-1] == '\n'):
        tokens[len(tokens) - 1] = tokens[len(tokens) - 1][0:-1]

    qcont = tokens[0]
    rcont = tokens[2]
    if (rcont == '*'):
        return -1, -1, -1, -1, "", ""

    l = int(tokens[3])
    cigar = tokens[5]

    cntSH = 0
    while (cntSH < len(cigar) and (cigar[cntSH].isdigit() or cigar[cntSH] == 'S' or cigar[cntSH] == 'H') ):
        cntSH += 1

    letCnt = parse_cigar(cigar[:cntSH])
    lq = letCnt['S'] + letCnt['H']

    letCnt = parse_cigar(cigar)

    rq = lq + letCnt['M'] + letCnt['='] + letCnt['X'] + letCnt['I']
    r = l + letCnt['M'] + letCnt['='] + letCnt['X'] + letCnt['D'] + letCnt['N']

    if ((int(tokens[1]) & (1 << 4)) != 0):
        rq, lq = lq, rq

    return lq, rq, l, r, qcont, rcont


def save_scaffolds_from_fasta(contig_file_name, lib):
    prevdir = os.getcwd()
    lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
    os.chdir(lib_dir)
    os.system("minimap2 -ax asm5 " + lib.path[0] + " " + contig_file_name + " > out.sam")
    #os.system("nucmer " + lib.path[0] + " " + contig_file_name)
    #os.system("show-coords out.delta -THrgl > out.coords")

    global cntedge
    global cntlib
    global idbyname

    output_json['libs'].append({'id': cntlib, 'color': lib.color, 'name': lib.label, 'type': 'SCAFF', 'scaffolds': []})

    contigsAlignment = dict()
    rcontlist = []

    scafflen_by_name = dict()
    save_lens_from_sam(scafflen_by_name, "out.sam")

    with open("out.sam") as g:
        for line in g:
            lq, rq, l, r, qcont, rcont = get_align_from_sam_line(line)
            if (lq == -1):
                continue

            chrlen = scafflen_by_name[rcont]
            if (lq > rq):
                qcont += "-rev"
                lq, rq = rq, lq

            id = idbyname[qcont]
            if (lq < 2 and rq > lenbyid[id] - 2):
                if (rcont not in contigsAlignment):
                    rcontlist.append(rcont)
                    rcontlist.append(rcont + "-rev")
                    contigsAlignment[rcont] = []
                    contigsAlignment[rcont + "-rev"] = []

                contigsAlignment[rcont].append((l, r, id))
                contigsAlignment[rcont + "-rev"].append((chrlen - r, chrlen - l, id^1))


    scafnum = 0
    for rc in rcontlist:
        contigsAlignment[rc].sort(key=lambda x: (x[0], -x[1]))
        output_json['libs'][-1]['scaffolds'].append({'name': rc, 'edges': []})

        lst = 0
        for i in range(1, len(contigsAlignment[rc])):
            if (contigsAlignment[rc][i][0] >= contigsAlignment[rc][lst][1] - 100):
                output_json['libs'][-1]['scaffolds'][-1]['edges'].append({})
                output_json['libs'][-1]['scaffolds'][-1]['edges'].append({'id': cntedge, 'from': contigsAlignment[rc][lst][2],
                                                                          'to': contigsAlignment[rc][i][2], 'weight': 1,
                                                                          'len': contigsAlignment[rc][i][0] - contigsAlignment[rc][lst][1]})
                cntedge += 1
                lst = i
        scafnum += 1

    cntlib += 1

    os.chdir(prevdir)


def save_scaffolds_from_gfa(lib):
    global cntedge
    global cntlib
    global idbyname
    with open(lib.path[0]) as g:
        output_json['libs'].append({'id': cntlib, 'color': lib.color, 'name': lib.label, 'type': 'SCAFF', 'scaffolds': []})

        scafnum = 0

        for line in g:
            tokens = line.split()
            if (tokens[0] != 'P'):
                continue

            output_json['libs'][-1]['scaffolds'].append({'name': tokens[1], 'edges': []})
            output_json['libs'][-1]['scaffolds'].append({'name': tokens[1] + "-rev", 'edges': []})

            nodeslist = []
            tt = tokens[2].split(',')
            for i in range(0, len(tt)):
                nm = tt[i][:-1]
                if (tt[i][-1] == '+'):
                    nodeslist.append(idbyname[nm])
                else:
                    nodeslist.append(idbyname[nm]^1)


            for i in range(1, len(nodeslist)):
                output_json['libs'][-1]['scaffolds'][-2]['edges'].append({'id': cntedge, 'from': nodeslist[i - 1],
                                                                          'to': nodeslist[i], 'weight': 1})
                cntedge += 1

            for i in range(len(nodeslist)-2, -1, -1):
                output_json['libs'][-1]['scaffolds'][-1]['edges'].append({'id': cntedge, 'from': nodeslist[i + 1]^1,
                                                                          'to': nodeslist[i]^1, 'weight': 1})
                cntedge += 1

            scafnum += 2
    cntlib += 1


def save_scaffolds_from_path(lib):
    global cntedge
    global cntlib
    global idbyname
    global partnametoname
    with open(lib.path[0]) as g:
        output_json['libs'].append({'id': cntlib, 'color': lib.color, 'name': lib.label, 'type': 'SCAFF', 'scaffolds': []})

        scafnum = 0
        lines = g.readlines()

        for il in range(0, len(lines), 2):
            scafname = lines[il]
            if (scafname[-1] == '\n'):
                scafname = scafname[:-1]
            if (scafname[-1] == "'"):
                continue

            tokens = lines[il + 1].split(",")
            if (tokens[len(tokens) - 1] == '\n'):
                tokens.pop()

            if (tokens[len(tokens) - 1][-1] == '\n'):
                tokens[len(tokens) - 1] = tokens[len(tokens) - 1][:-1]

            output_json['libs'][-1]['scaffolds'].append({'name': scafname, 'edges': []})
            output_json['libs'][-1]['scaffolds'].append({'name': scafname + "-rev", 'edges': []})

            nodeslist = []
            for i in range(0, len(tokens)):
                nm = partnametoname[tokens[i][:-1]]
                if (tokens[i][-1] == '+'):
                    nodeslist.append(idbyname[nm])
                else:
                    nodeslist.append(idbyname[nm]^1)


            for i in range(1, len(nodeslist)):
                output_json['libs'][-1]['scaffolds'][-2]['edges'].append({'id': cntedge, 'from': nodeslist[i - 1],
                                                                          'to': nodeslist[i], 'weight': 1})
                cntedge += 1

            for i in range(len(nodeslist)-2, -1, -1):
                output_json['libs'][-1]['scaffolds'][-1]['edges'].append({'id': cntedge, 'from': nodeslist[i + 1] ^ 1,
                                                                          'to': nodeslist[i] ^ 1, 'weight': 1})
                cntedge += 1

            scafnum += 2
    cntlib += 1


def save_scaffolds(contig_file_name, args):
    for lib in args.libs["scafinfo"]:
        save_scaffolds_from_info(lib)

    for lib in args.libs["scafpath"]:
        save_scaffolds_from_path(lib)

    for lib in args.libs["scaffolds"]:
        save_scaffolds_from_fasta(contig_file_name, lib)

    for lib in args.libs["gfa"]:
        save_scaffolds_from_gfa(lib)


def get_lib_pos_by_id(libid):
    for i in range(len(output_json['libs'])):
        if (output_json['libs'][i]['id'] == libid):
            return i
    return -1

def add_conection_to_res_file():
    global cntlib
    global cntedge
    if (not os.path.isfile("graph.gr")):
        return

    with open("graph.gr") as g:
        cntlib = int(g.readline())

        for i in range(cntlib):
            libsinfo = g.readline().split(" ")
            libsinfo[4] = libsinfo[4][:-1]
            output_json['libs'].append({'id': int(libsinfo[1]), 'color': libsinfo[2], 'name': libsinfo[3], 'type': libsinfo[4], 'edges': []})

        nodecnt = int(g.readline())

        for i in range(nodecnt):
            g.readline()

        cntedge = int(g.readline())
        for i in range(cntedge):
            curs = g.readline()
            extraInfo = ""
            if ("\""  in curs):
                extraInfo = curs.split("\"")[1]
            edgesinfo = curs.split()

            output_json['libs'][get_lib_pos_by_id(int(edgesinfo[4]))]['edges'].append({'id': int(edgesinfo[1]), 'from': int(edgesinfo[2]),
                                                                                  'to': int(edgesinfo[3]), 'weight': float(edgesinfo[5]),
                                                                                  'len': edgesinfo[6], 'info': extraInfo})


def add_refcoord_to_res_file():
    if (len(args.refcoord) == 0):
        return

    lib = args.refcoord[0]

    chrid = {}
    chrlen = []
    fasta_seq = SeqIO.parse(open(lib[0]), 'fasta')
    curid = 0

    for fasta in fasta_seq:
        name, lenn = fasta.id, len(fasta.seq.tostring())
        chrid[name] = curid
        chrid[name + "-rev"] = curid + 1
        chrlen.append(lenn)
        chrlen.append(lenn)
        output_json["chromosomes"].append({"id": curid, "name": name, "len": lenn})
        output_json["chromosomes"].append({"id": curid + 1, "name": name + "-rev", "len": lenn})
        curid += 2

    global idbyname
    global lenbyid

    lastname = '-'
    #TODO: del file
    g = open("out.coords", "w")

    with open(lib[1]) as cf:
        for line in cf:
            if ("[S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]" in line or "====" in line):
                continue
            info = line.split(" ")
            info[-1] = info[-1][:-1]
            vid = idbyname[info[12]]
            curid = chrid[info[11].split('_')[0]]
            lq = int(info[3])
            rq = int(info[4])
            l = int(info[0])
            r = int(info[1])
            lenf = chrlen[curid]
            if ((max(rq, lq) - min(rq, lq)) * 100 < lenbyid[vid]):
                continue
            if (lq > rq):
                vid ^= 1

            output_json["alignments"].append({"coord_begin": l, "coord_end": r,
                                              "chr_id": curid, "node_id": vid})
            output_json["alignments"].append({"coord_begin": lenf - r, "coord_end": lenf - l,
                                              "chr_id": curid + 1, "node_id": vid ^ 1})
            g.write(str(l) + " " + str(r) + " " + str(lq) + " " + str(rq) + " 0 0 0 " + str(lenf) + " 0 " + info[11].split('_')[0] + " " + info[12] + "\n")
    g.close()

    
def getRefFileName(fileName):
    return fileName.split('/')[-1].split('.')[0]


def merge_ref_files(ref_libs):
    with open("ref_merge.fasta", "w") as out:
        for i in range(len(ref_libs)):
            for record in SeqIO.parse(ref_libs[i], "fasta"):
                record.id = getRefFileName(ref_libs[i]) + "_" + record.id
                SeqIO.write(record, out, "fasta")

    return Lib([os.path.abspath("ref_merge.fasta")], "ref", "ref")

def add_ref_to_res_file(contig_file_name):
    if (len(args.ref) == 0):
        return

    lib = merge_ref_files(args.ref)

    prevdir = os.getcwd()
    lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
    os.chdir(lib_dir)
    os.system("minimap2 -ax asm5 " + lib.path[0] + " " + contig_file_name + " > out.sam")
    #os.system("nucmer -b 10000 " + lib.path[0] + " " + contig_file_name)
    #os.system("show-coords out.delta -THrgl > out.coords")

    global idbyname
    global lenbyid

    chrnameToId = dict()

    curid = -2

    chrm_len_by_name = dict()
    save_lens_from_sam(chrm_len_by_name, "out.sam")

    with open("out.sam") as cf:
        for line in cf:
            lq, rq, l, r, qcont, chrname = get_align_from_sam_line(line)
            if (lq == -1):
                continue
            vid = idbyname[qcont]
            lenf = int(chrm_len_by_name[chrname])
            if (chrname not in chrnameToId):
                curid += 2
                output_json["chromosomes"].append({"id": curid, "name": chrname, "len": lenf})
                output_json["chromosomes"].append({"id": curid + 1, "name": chrname + "-rev", "len": lenf})
                chrnameToId[chrname] = curid

            if ((max(rq, lq) - min(rq, lq)) * 100 < lenbyid[vid]):
                continue
            
            if (lq > rq):
                vid ^= 1

            output_json["alignments"].append({"coord_begin": l, "coord_end": r, "chr_id": chrnameToId[chrname], "node_id": vid})
            output_json["alignments"].append({"coord_begin": lenf - r, "coord_end": lenf - l, "chr_id": chrnameToId[chrname] + 1, "node_id": vid^1})

    os.chdir(prevdir)


def merge_contigs(contigs):
    with open("contigs_merge.fasta", "w") as out:
        for i in range(len(contigs)):
            for record in SeqIO.parse(contigs[i], "fasta"):
                SeqIO.write(record, out, "fasta")

    return os.path.abspath("contigs_merge.fasta")

def fastg_to_contigs(args):
    for lib in args.libs['fastg']:
        lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
        if not os.path.exists(lib_dir):
            os.makedirs(lib_dir)
        prevdir = os.getcwd()
        os.chdir(lib_dir)

        with open("contigs.fasta", "w") as out:
            for record in SeqIO.parse(lib.path[0], "fasta"):
                record.id = record.id.split(':')[0].split(';')[0]
                if (record.id[-1] != '\''):
                    partnametoname[record.id.split('_')[1]]=record.id
                    SeqIO.write(record, out, "fasta")

        if (args.contigs == None):
            args.contigs = []

        args.contigs.append(os.path.abspath('contigs.fasta'))
        os.chdir(prevdir)


def gfa_to_contigs(args):
    for lib in args.libs['gfa']:
        lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
        if not os.path.exists(lib_dir):
            os.makedirs(lib_dir)
        prevdir = os.getcwd()
        os.chdir(lib_dir)

        with open("contigs.fasta", "w") as out:
            lines = [line.rstrip('\n') for line in open(lib.path[0])]
            for line in lines:
                parts = line.split()
                if (parts[0] == 'S'):
                    record = SeqRecord(Seq(parts[2], IUPAC.ambiguous_dna), id=parts[1], description='')
                    SeqIO.write(record, out, "fasta")

        if (args.contigs == None):
            args.contigs = []

        args.contigs.append(os.path.abspath('contigs.fasta'))
        os.chdir(prevdir)


def run(args):
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    if args.contigs == None and (('libs' not in args) or (len(args.libs['fastg']) == 0)) and (('libs' not in args) or (len(args.libs['gfa']) == 0)):
        log.err("none contig/FASTG/GFA file provide")
        sys.exit()

    if (args.contigs != None):
        for i in range(len(args.contigs)):
            args.contigs[i] = os.path.abspath(args.contigs[i][0])

    main_out_dir = os.path.abspath(".") + "/"

    if args.local_output_dir != None:
        main_out_dir = os.path.abspath(args.local_output_dir[0]) + "/"

    if args.graph != None:
        for i in range(len(args.graph)):
            args.graph[i] = os.path.abspath(args.graph[i][0])

    if args.ref != None:
        for i in range(len(args.ref)):
            args.ref[i] = os.path.abspath(args.ref[i][0])

    if args.refcoord != None:
        for i in range(len(args.refcoord)):
            args.refcoord[i][0] = os.path.abspath(args.refcoord[i][0])
            args.refcoord[i][1] = os.path.abspath(args.refcoord[i][1])


    out_dir = main_out_dir + "tmp/"
    log.log("OUTPUT DIR: " + out_dir)
    directory = os.path.dirname(out_dir)
    if not os.path.exists(directory):
        log.log("MKDIR")
        os.makedirs(directory)
    os.chdir(directory)

    if 'libs' in args:
        fastg_to_contigs(args)
        gfa_to_contigs(args)

    contig_file_name = ""
    if (len(args.contigs) == 1):
        contig_file_name = args.contigs[0]
    else:
        contig_file_name = merge_contigs(args.contigs)

    if args.color != None and len(args.color) != args.lib_cnt:
        log.err("wrong number of color provide, lib cnt = " + str(args.lib_cnt) + " color cnt = " + str(len(args.color)))

    if args.label != None and len(args.label) != args.lib_cnt:
        log.err("wrong number of labels provide")

    if 'libs' in args:
        for lib_type in libsType:
            for lib in args.libs[lib_type]:
                lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
                if not os.path.exists(lib_dir):
                    os.makedirs(lib_dir)

                if args.color != None:
                    lib.color = args.color[lib.id]
                if args.label != None:
                    lib.label = args.label[lib.id]

        if (not os.path.exists(os.path.dirname(os.path.abspath("ref_ref") + "/"))):
            os.makedirs(os.path.dirname(os.path.abspath("ref_ref") + "/"))
        alig_reads(contig_file_name, args)
        build_graph(contig_file_name, args)

    merge_graph(args) #result file graph.gr

    gen_id_from_contig_file(contig_file_name)
    output_json['libs'] = []
    output_json['chromosomes'] = []
    output_json["alignments"] = []

    add_conection_to_res_file()
    if 'libs' in args:
        save_scaffolds(contig_file_name, args)

    if args.ref != None:
        add_ref_to_res_file(contig_file_name)

    if args.refcoord != None:
        add_refcoord_to_res_file()

    f = open("data.json", 'w')
    f.write(json.dumps(output_json))
    f.close()

    directory = os.path.dirname(main_out_dir)
    os.chdir(directory)

    os.system("cp -r " + path_to_exec_dir + "/scripts ./")
    move("tmp/data.json", "./scripts/data.json")
    os.system("cp " + path_to_exec_dir + "/mainPage.html ./main.html")
    return

args = parse_args()
run(args)
