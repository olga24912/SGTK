#!/usr/bin/python3
import sys
import os
import argparse
from shutil import copyfile
from shutil import move
from Bio import SeqIO

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

libsType = {"rnap", "rnas", "dnap", "ref", "scafinfo", "scaffolds"}

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
    parser.add_argument("--scaffolds", "-s", nargs=1, dest="scaffolds", help="path to scaffolds in fasta format", type=str, action=StoreArgAction)
    parser.add_argument("--rna-p", dest="rnap", nargs=2, help="path to rna pair reads file", type=str, action=StoreArgAction)
    parser.add_argument("--rna-s", dest="rnas", nargs=1, help="path to rna read file", type=str, action=StoreArgAction)
    parser.add_argument("--local_output_dir", "-o", nargs=1, help="use this output dir", type=str)
    parser.add_argument("--ref", dest="ref", nargs=1, help="path to reference", type=str, action=StoreArgAction)
    parser.add_argument("--dna-p", dest="dnap", nargs=2, help="path to dna pair reads file", type=str, action=StoreArgAction)
    parser.add_argument("--scafinfo", nargs=1, help="path to .info file with info about scaffolds", type=str, action=StoreArgAction)
    parser.add_argument("--label", "-l", nargs='*', help="list with labels for all libs in definition order", type=str, action='store')
    parser.add_argument("--color", nargs='*', help="list with color for all libs in definition order", type=str, action='store')
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

    os.system("STAR --runThreadN 20 --genomeDir ../genomeDir --outReadsUnmapped Fastx --readFilesIn " + rnap.path[0])
    os.system("mv Aligned.out.sam rna1.sam")
    if rnap.path[0][-1] == "q":
        os.system("mv Unmapped.out.mate1 Unmapped1.fastq")
        unm1 = "../" + rnap.name + "/Unmapped1.fastq"
    else:
        os.system("mv Unmapped.out.mate1 Unmapped1.fasta")
        unm1 = "../" + rnap.name + "/Unmapped1.fasta"

    os.system("STAR --runThreadN 20 --genomeDir ../genomeDir --outReadsUnmapped Fastx --readFilesIn " + rnap.path[1])
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

    os.system("bowtie2-build " + contig_file_name + " contig")
    if dnap.path[0].endswith(".fa") or dnap.path[0].endswith(".fasta") or dnap.path[0].endswith(".mfa") or dnap.path[0].endswith(".fna"):
        os.system("bowtie2 -x contig -f --ignore-quals -U " + dnap.path[0] + " -S dna1.sam")
        os.system("bowtie2 -x contig -f --ignore-quals -U " + dnap.path[1] + " -S dna2.sam")
    else:
        os.system("bowtie2 -x contig -U " + dnap.path[0] + " -S dna1.sam")
        os.system("bowtie2 -x contig -U " + dnap.path[1] + " -S dna2.sam")
    os.chdir(prevdir)

def alig_single_rna_reads(rnas):
    prevdir = os.getcwd()
    lib_name = rnas.name + "_50"
    lib_dir = os.path.dirname(os.path.abspath(lib_name) + "/")
    if not os.path.exists(lib_dir):
        os.makedirs(lib_dir)
    os.chdir(lib_name)
    unm = ""
    os.system("STAR --runThreadN 20 --genomeDir ../genomeDir --outReadsUnmapped Fastx --readFilesIn " + rnas)
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
        os.system("STAR --runMode genomeGenerate --genomeDir genomeDir --runThreadN 20 --genomeSAindexNbases 10 --genomeFastaFiles " +
                  contig_file_name + " --limitGenomeGenerateRAM 90000000000")
    except:
        log.err(sys.exc_info()[0])
        return

    for i in range(len(args.libs["rnap"])):
        alig_pair_rna_reads(args.libs["rnap"][i])

    for i in range(len(args.libs["dnap"])):
        alig_pair_dna_reads(args.libs["dnap"][i], contig_file_name)

    for i in range(len(args.libs["rnas"])):
        alig_single_rna_reads(args.libs["rnas"][i])

    return

def runGraphBuilder(lib_name, prevdir, type, label):
    log.log("START BUILD GRAPH: " + lib_name)
    lib_dir = os.path.dirname(os.path.abspath(lib_name) + "/")
    os.chdir(lib_dir)
    os.system(path_to_exec_dir + "build " + type + " rna1.sam rna2.sam " + label)
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

    for lib in args.libs["dnap"]:
        prevdir = os.getcwd()
        log.log("START BUILD GRAPH: " + lib.label)
        lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
        os.chdir(lib_dir)
        os.system(path_to_exec_dir + "build DNA_PAIR dna1.sam dna2.sam 1000000000 " + lib.label)
        os.chdir(prevdir)

    for lib in args.libs["ref"]:
        prevdir = os.getcwd()
        log.log("START BUILD GRAPH: " + lib.label)
        lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
        os.chdir(lib_dir)
        os.system(path_to_exec_dir + "build REF " + lib.path[0] + " " + contig_file_name + " 500 " + lib.label)
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
        os.system(path_to_exec_dir + "filter " + os.path.abspath("../filter_config"))
        os.chdir(prevdir)
    return

def merge_graph(args):
    merge_lib_rna(args.libs["rnap"])

    merge_list = ""

    for lib_type in libsType:
        for lib in args.libs[lib_type]:
            if lib_type != "scafinfo" and lib_type != "scafpath":
                lib.fix_graph_file()
                if lib_type != "rnas":
                    merge_list += lib.name + "/graph.gr "
                else:
                    merge_list += lib.name + "_50/graph.gr"
                    merge_list += lib.name + "_30/graph.gr"

    merge_list += "graph.gr"
    os.system(path_to_exec_dir + "mergeGraph " + merge_list)
    return

idbyname = dict()
lenbyid = []
cntedge = 0
cntlib = 0

#save id to idbyname
#write node info to data
def gen_id_from_contig_file(contig_file_name, f):
    fasta_seq = SeqIO.parse(open(contig_file_name), 'fasta')
    id = 0
    nodestr = "var scaffoldnodes = ["
    for fasta in fasta_seq:
        name, lenn = fasta.id, len(fasta.seq.tostring())
        idbyname[name] = id
        idbyname[name + "-rev"] = id + 1
        lenbyid.append(lenn)
        lenbyid.append(lenn)
        if (id != 0):
            nodestr += ', '
        nodestr += "new ScaffoldNode(" + str(id) + ", '" + name + "', " + str(lenn) + "), "
        nodestr += "new ScaffoldNode(" + str(id + 1) + ", '" + name + "-rev', " + str(lenn) + ")"
        id += 2
    nodestr += "];"
    f.write(nodestr)


def save_scaffolds_from_info(lib, f):
    global cntedge
    global cntlib
    global idbyname
    with open(lib.path[0]) as g:
        f.write("scaffoldlibs.push(new ScaffoldEdgeLib(" + str(cntlib) + ", '" + str(lib.color) + "', '" + str(lib.label) + "', 'SCAFF'));\n")

        scafnum = 0

        for line in g:
            tokens = line.split(" ")
            if (tokens[len(tokens) - 1] == '\n'):
                tokens.pop()

            f.write("scaffoldlibs["+ str(cntlib) +"].scaffolds.push(new Scaffold('" + tokens[0][1:] + "'));\n")
            f.write("scaffoldlibs["+ str(cntlib) +"].scaffolds.push(new Scaffold('" + tokens[0][1:] + "-rev'));\n")

            nodeslist = []
            for i in range(1, len(tokens), 3):
                nm = tokens[i][1:]
                if (tokens[i + 2][0] == '+'):
                    nodeslist.append(idbyname[nm])
                else:
                    nodeslist.append(idbyname[nm]^1)


            for i in range(1, len(nodeslist)):
                f.write("scaffoldedges.push(new ScaffoldEdge(" + str(cntedge) + ", "+ str(nodeslist[i - 1]) +
                        ", " + str(nodeslist[i]) + ", " + str(cntlib) + ", 1));\n")
                f.write("scaffoldedges["+str(cntedge)+"].name='"+ tokens[0][1:] + "';\n")
                f.write("scaffoldlibs["+ str(cntlib) +"].scaffolds["+str(scafnum) +"].edges.push(scaffoldedges["+str(cntedge)+"]);\n")
                cntedge += 1

            for i in range(len(nodeslist)-2, -1, -1):
                f.write("scaffoldedges.push(new ScaffoldEdge(" + str(cntedge) + ", "+ str(nodeslist[i + 1]^1) +
                        ", " + str(nodeslist[i]^1) + ", " + str(cntlib) + ", 1));\n")
                f.write("scaffoldedges["+str(cntedge)+"].name='"+ tokens[0][1:] + "-rev';\n")
                f.write("scaffoldlibs["+ str(cntlib) +"].scaffolds["+str(scafnum+1) +"].edges.push(scaffoldedges["+str(cntedge)+"]);\n")
                cntedge += 1

            scafnum += 2
    cntlib += 1


def save_scaffolds_from_fasta(contig_file_name, lib, f):
    prevdir = os.getcwd()
    lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
    os.chdir(lib_dir)
    os.system("nucmer " + lib.path[0] + " " + contig_file_name)
    os.system("show-coords out.delta -THrgl > out.coords")

    global cntedge
    global cntlib
    global idbyname
    with open("out.coords") as g:
        f.write("scaffoldlibs.push(new ScaffoldEdgeLib(" + str(cntlib) + ", '" + str(lib.color) + "', '" + str(lib.label) + "', 'SCAFF'));\n")

        contigsAlignment = dict()
        rcontlist = []

        for line in g:
            tokens = line.split("\t")
            print(tokens)
            if (tokens[len(tokens) - 1] == '\n'):
                tokens.pop()
            if (tokens[len(tokens) - 1][-1] == '\n'):
                tokens[len(tokens) - 1] = tokens[len(tokens) - 1][0:-1]

            lq = int(tokens[2])
            rq = int(tokens[3])
            l = int(tokens[0])
            r = int(tokens[1])
            qcont = tokens[10]
            rcont = tokens[9]
            chrlen = int(tokens[7])
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


        print(rcontlist)
        scafnum = 0
        for rc in rcontlist:
            contigsAlignment[rc].sort()
            print(rc)
            print(contigsAlignment[rc])
            f.write("scaffoldlibs["+ str(cntlib) +"].scaffolds.push(new Scaffold('" + rc + "'));\n")

            lst = 0
            for i in range(1, len(contigsAlignment[rc])):
                if (contigsAlignment[rc][i][0] >= contigsAlignment[rc][lst][1]):
                    f.write("scaffoldedges.push(new ScaffoldEdge(" + str(cntedge) + ", "+ str(contigsAlignment[rc][lst][2]) +
                        ", " + str(contigsAlignment[rc][i][2]) + ", " + str(cntlib) + ", 1));\n")
                    f.write("scaffoldedges["+str(cntedge)+"].name='"+ rc + "';\n")
                    f.write("scaffoldlibs["+ str(cntlib) +"].scaffolds["+str(scafnum) +"].edges.push(scaffoldedges["+str(cntedge)+"]);\n")
                    cntedge += 1
                    lst = i

            scafnum += 1

        cntlib += 1

    os.chdir(prevdir)

def save_scaffolds(contig_file_name, args, f):
    for lib in args.libs["scafinfo"]:
        save_scaffolds_from_info(lib, f)

    for lib in args.libs["scaffolds"]:
        save_scaffolds_from_fasta(contig_file_name, lib, f)


def add_conection_to_res_file(f):
    with open("graph.gr") as g:
        cntlib = int(g.readline())

        for i in range(cntlib):
            libsinfo = g.readline().split(" ")
            libsinfo[4] = libsinfo[4][:-1]
            f.write("scaffoldlibs.push(new ScaffoldEdgeLib(" + libsinfo[1] + ", '" + libsinfo[2] + "', '" + libsinfo[3] + "', '" + libsinfo[4] + "'));\n")

        nodecnt = int(g.readline())

        for i in range(nodecnt):
            g.readline()

        edgescnt = int(g.readline())
        for i in range(edgescnt):
            edgesinfo = (g.readline()).split(" ")
            f.write("scaffoldedges.push(new ScaffoldEdge(" + edgesinfo[1] + ", " + edgesinfo[2] + ", " + edgesinfo[3] + ", " + edgesinfo[4] + ", " + edgesinfo[5] + "));\n")


def add_ref_to_res_file(contig_file_name, f):
    lib = args.libs["ref"][0]

    prevdir = os.getcwd()
    lib_dir = os.path.dirname(os.path.abspath(lib.name) + "/")
    os.chdir(lib_dir)
    os.system("nucmer " + lib.path[0] + " " + contig_file_name)
    os.system("show-coords out.delta -THrgl > out.coords")

    global idbyname

    chrlist = []
    chralig = []

    curid = -2
    lastname = '-'

    with open("out.coords") as cf:
        for line in cf:
            info = line.split("\t")
            info[10] = info[10][:-1]
            vid = idbyname[info[10]]
            chrname = info[9]
            lenf = int(info[7])
            if (chrname != lastname):
                curid += 2
                chrlist.append("new Chromosome(" + str(curid) + ", '" + chrname + "', " + str(lenf) + ")")
                chrlist.append("new Chromosome(" + str(curid + 1) + ", '" + chrname + "-rev', " + str(lenf) + ")")
                chralig.append([])
                chralig.append([])
                lastname = chrname

            lq = int(info[2])
            rq = int(info[3])
            l = int(info[0])
            r = int(info[1])

            if (lq > rq):
                vid ^= 1

            chralig[curid].append("new Alignment(" + str(l) + ", " + str(r) + ", " + str(curid) + ", " + str(vid) + ")")
            chralig[curid + 1].append("new Alignment(" + str(lenf - r) + ", " + str(lenf - l) + ", " + str(curid + 1) + ", " + str(vid^1) + ")")

    f.write("var chromosomes = [")
    for i in range(len(chrlist)):
        f.write(chrlist[i])
        if (i != len(chrlist) - 1):
            f.write(", ")
    f.write("];\n")

    for i in range(len(chrlist)):
        f.write("chromosomes[" + str(i) + "].alignments = " + "[")
        for j in range(len(chralig[i])):
            f.write(chralig[i][j])
            if (j != len(chralig[i]) - 1):
                f.write(", ")
        f.write("];\n")

    f.close()

    os.chdir(prevdir)


def run(args):
    if args.contigs == None:
        log.err("none contig file provide")
        return

    contig_file_name = os.path.abspath(args.contigs[0])

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
    merge_graph(args) #result file graph.gr

    f = open("data.js", 'w')
    gen_id_from_contig_file(contig_file_name, f)
    f.write("var scaffoldlibs = [];\n")
    f.write("var scaffoldedges = [];\n")
    add_conection_to_res_file(f)
    save_scaffolds(contig_file_name, args, f)
    add_ref_to_res_file()
    f.write("var scaffoldgraph = new ScaffoldGraph(scaffoldlibs, scaffoldnodes, scaffoldedges);\n")
    f.close()

    directory = os.path.dirname(main_out_dir)
    os.chdir(directory)

    move("tmp/data.js", "data.js")
    return

args = parse_args()
run(args)