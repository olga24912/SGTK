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
    parser.add_argument("--local_output_dir", "-o", nargs=1, help="output folder", type=str)
    parser.add_argument("--contigs", "-c", nargs=1, dest="contigs", help="path to contigs  in FASTA format", type=str)
    parser.add_argument("--rna-p", dest="rnap", nargs=2, help="paths to RNA-Seq paired reads  in SAM/BAM format", type=str, action='append')
    parser.add_argument("--rna-s", dest="rnas", nargs=1, help="path to RNA-Seq single reads  in SAM/BAM format", type=str,action='append')
    parser.add_argument("--gene_annotation", nargs=1, help="path to gff file with gene anotation for inputs contigs (optional)", type=str)

    args = parser.parse_args()
    return args


def alig_split(lib_name, reads, flag, checkpoints, cpf):
    if not (lib_name + " alig" in checkpoints):
        prevdir = os.getcwd()
        log.log("align " + lib_name)
        lib_dir = os.path.dirname(os.path.abspath(lib_name) + "/")
        if not os.path.exists(lib_dir):
            os.makedirs(lib_dir)
        os.chdir(lib_dir)

        os.system(path_to_exec_dir + "readSplitter " + str(flag) + " " + reads + " reads1.fasta reads2.fasta")
        os.system("STAR --runThreadN 20 --genomeDir ../genomeDir --readFilesIn reads1.fasta --outSAMtype BAM Unsorted")
        os.system("samtools sort -n Aligned.out.bam -o rna1.bam")
        os.system("rm -r _STARtmp")
        os.system("STAR --runThreadN 20 --genomeDir ../genomeDir --readFilesIn reads2.fasta --outSAMtype BAM Unsorted")
        os.system("samtools sort -n Aligned.out.bam -o rna2.bam")
        os.system("rm -r _STARtmp")

        os.chdir(prevdir)
    else:
        log.warn("skip step align " + lib_name +
                 ". If you don't want skip this step delete raw \"" + lib_name + " alig" "\" from checkpoints file")

    cpf.write(lib_name + " alig\n")


def alig_pair_reads(i, rnap1, rnap2, checkpoints, cpf):
    prevdir = os.getcwd()
    lib_name = "rnap" + str(i)
    unm1 = ""
    unm2 = ""
    if not(lib_name + " alig" in checkpoints):
        log.log("align " + lib_name + ": " + rnap1 + " " + rnap2)
        lib_dir = os.path.dirname(os.path.abspath(lib_name) + "/")
        if not os.path.exists(lib_dir):
            os.makedirs(lib_dir)
        os.chdir(lib_dir)

        os.system("STAR --runThreadN 20 --genomeDir ../genomeDir --outReadsUnmapped Fastx --readFilesIn " + rnap1 +
                  " --outSAMtype BAM Unsorted")
        os.system("samtools sort -n Aligned.out.bam -o rna1.bam")
        os.system("rm -r _STARtmp")
        if rnap1[-1] == "q":
            os.system("mv Unmapped.out.mate1 Unmapped1.fastq")
            unm1 = "../" + lib_name + "/Unmapped1.fastq"
        else:
            os.system("mv Unmapped.out.mate1 Unmapped1.fasta")
            unm1 = "../" + lib_name + "/Unmapped1.fasta"

        os.system("STAR --runThreadN 20 --genomeDir ../genomeDir --outReadsUnmapped Fastx --readFilesIn " + rnap2 +
                  " --outSAMtype BAM Unsorted")
        os.system("samtools sort -n Aligned.out.bam -o rna2.bam")
        os.system("rm -r _STARtmp")
        if rnap2[-1] == "q":
            os.system("mv Unmapped.out.mate1 Unmapped2.fastq")
            unm2 = "../" + lib_name + "/Unmapped2.fastq"
        else:
            os.system("mv Unmapped.out.mate1 Unmapped2.fasta")
            unm2 = "../" + lib_name + "/Unmapped2.fasta"

        os.chdir(prevdir)
    else:
        log.warn("skip step align reads: " + rnap1 + " and " + rnap2 +
                 ". If you don't want skip this step delete raw \"" + lib_name + " alig" "\" from checkpoints file")

        if (os.path.isfile("../" + lib_name + "/Unmapped1.fasta")):
            unm1 = "../" + lib_name + "/Unmapped1.fasta"
        else:
            unm1 = "../" + lib_name + "/Unmapped1.fastq"

        if (os.path.isfile("../" + lib_name + "/Unmapped2.fasta")):
            unm2 = "../" + lib_name + "/Unmapped2.fasta"
        else:
            unm2 = "../" + lib_name + "/Unmapped2.fastq"

    cpf.write(lib_name + " alig\n")
    alig_split("rnap" + str(i) + "_50_1", unm1, 0, checkpoints, cpf)
    alig_split("rnap" + str(i) + "_50_2", unm2, 0, checkpoints, cpf)
    alig_split("rnap" + str(i) + "_30_1", "../" + lib_name + "/rna1.bam", 1, checkpoints, cpf)
    alig_split("rnap" + str(i) + "_30_2", "../" + lib_name + "/rna2.bam", 1, checkpoints, cpf)


def alig_single_reads(i, rnas, checkpoints, cpf):
    prevdir = os.getcwd()
    lib_name = "rnas" + str(i) + "_50"
    if not(lib_name + " alig" in checkpoints):
        log.log("align " + rnas + " reads to detect unalign reads")
        lib_dir = os.path.dirname(os.path.abspath(lib_name) + "/")
        if not os.path.exists(lib_dir):
            os.makedirs(lib_dir)
        os.chdir(lib_name)
        unm = ""
        os.system("STAR --runThreadN 20 --genomeDir ../genomeDir --outReadsUnmapped Fastx --readFilesIn " + rnas +
                  " --outSAMtype BAM Unsorted")
        os.system("samtools sort -n Aligned.out.bam -o rna.bam")
        os.system("rm -r _STARtmp")
        if rnas[-1] == "q":
            os.system("mv Unmapped.out.mate1 Unmapped.fastq")
            unm = "Unmapped.fastq"
        else:
            os.system("mv Unmapped.out.mate1 Unmapped.fasta")
            unm = "Unmapped.fasta"


        os.chdir(prevdir)
        alig_split(lib_name, unm, 0, checkpoints, cpf)
    else:
        log.warn("skip step align " + lib_name +
                 ". If you don't want skip this step delete raw \"" + lib_name + " alig" "\" from checkpoints file")

    cpf.write(lib_name + " alig")
    alig_split("rnas" + str(i) + "_30", "../rnas" + str(i) + "_30/rna.bam", 1, checkpoints, cpf)


def alig_reads(contig_file_name, rnap_list, rnas_list, checkpoints, cpf):
    log.log("PHASE 1: READS' ALIGNMENT")

    genome_dir = "genomeDir"
    gen_dir = os.path.dirname(os.path.abspath(genome_dir) + "/")
    if not os.path.exists(gen_dir):
        os.makedirs(gen_dir)

    if not("geneDir" in checkpoints):
        log.log("genome generation for contigs: " + contig_file_name)
        try:
            os.system("STAR --runMode genomeGenerate --genomeDir genomeDir --runThreadN 20 --genomeSAindexNbases 10 --genomeFastaFiles " +
                  contig_file_name + " --limitGenomeGenerateRAM 90000000000")
        except:
            log.err(sys.exc_info()[0])
            return
    else:
        log.warn("skip genome generation step. If you have not genomeDir directory delete raw \"geneDir\" from checkpoints file")

    cpf.write("geneDir\n")

    for i in range(len(rnap_list)):
        alig_pair_reads(i, rnap_list[i][0], rnap_list[i][1], checkpoints, cpf)

    for i in range(len(rnas_list)):
        alig_single_reads(i, rnas_list[i][0], checkpoints, cpf)

    log.log("FINISH PHASE 1")
    return


def runGraphBuilder(lib_name, prevdir, type, checkpoints, cpf):
    if (not (lib_name + " build" in checkpoints)):
        lib_dir = os.path.dirname(os.path.abspath(lib_name) + "/")
        os.chdir(lib_dir)
        os.system(path_to_exec_dir + "build " + type + " rna1.bam rna2.bam " + lib_name)
        os.chdir(prevdir)
    else:
        log.warn("skip building graph for " + lib_name +
                 ". If you don't want skip this step delete raw \"" + lib_name + " build" "\" from checkpoints file")
    cpf.write(lib_name + " build\n")
    return


def build_graph(contig_file_name, rnap_list, rnas_list, checkpoints, cpf):
    log.log("PHASE 2: GRAPH BUILDING")

    for i in range(len(rnap_list)):
        prevdir = os.getcwd()
        runGraphBuilder("rnap" + str(i), prevdir, "RNA_PAIR", checkpoints, cpf)
        runGraphBuilder("rnap" + str(i) + "_50_1", prevdir, "RNA_SPLIT_50", checkpoints, cpf)
        runGraphBuilder("rnap" + str(i) + "_50_2", prevdir, "RNA_SPLIT_50", checkpoints, cpf)
        runGraphBuilder("rnap" + str(i) + "_30_1", prevdir, "RNA_SPLIT_30", checkpoints, cpf)
        runGraphBuilder("rnap" + str(i) + "_30_2", prevdir, "RNA_SPLIT_30", checkpoints, cpf)

    for i in range(len(rnas_list)):
        prevdir = os.getcwd()
        runGraphBuilder("rnas" + str(i) + "_50", prevdir, "RNA_SPLIT_50", checkpoints, cpf)
        runGraphBuilder("rnas" + str(i) + "_30", prevdir, "RNA_SPLIT_30", checkpoints, cpf)

    log.log("FINISH PHASE 2")
    return


def merge_graph(rnap_list, rnas_list):
    log.log("PHASE 3: GRAPHS' MERGING")

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
    log.log("FINISH PHASE 3")
    return

def create_scaffolds(contig_file_name, rnap_list, rnas_list, exon_block_file_name):
    log.log("PHASE 4: GRAPH SIMPLIFICATION and SCAFFOLDS' CONSTRACTION")

    f = open("filter_config", 'w')
    f.write("uploadGraph out.gr\n")
    f.write("minContig 500\n")
    f.write("setExonBlock " + exon_block_file_name + "\n")
    f.write("mergeSimplePath " + contig_file_name + " scaffolds.fa\n")
    f.write("exit\n")
    f.close()

    os.system(path_to_exec_dir + "filter " + os.path.abspath("filter_config"))

    log.log("FINISH PHASE 4")
    return

def run(args):
    if args.contigs == None:
        log.err("none contig file provide")
        return

    contig_file_name = os.path.abspath(args.contigs[0])
    exon_block_file_name = ""
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

    if args.gene_annotation != None:
        exon_block_file_name = os.path.abspath(args.gene_annotation[0])
    else:
        log.err("none gene annotation file provide")
        return

    out_dir = main_out_dir + "tmp/"
    directory = os.path.dirname(out_dir)
    if not os.path.exists(directory):
        os.makedirs(directory)
    os.chdir(directory)

    filename = "checkpoints"
    checkpoints = [""]
    if os.path.isfile(filename):
        with open(filename) as f:
            checkpoints = f.read().splitlines()

    cpf = open('checkpoints', 'w')


    alig_reads(contig_file_name, rnap_list, rnas_list, checkpoints, cpf)
    build_graph(contig_file_name, rnap_list, rnas_list, checkpoints, cpf)
    merge_graph(rnap_list, rnas_list)
    create_scaffolds(contig_file_name, rnap_list, rnas_list, exon_block_file_name)

    directory = os.path.dirname(main_out_dir)
    os.chdir(directory)

    cpf.close()

    copyfile("tmp/scaffolds.fa", "scaffolds.fa")
    copyfile("tmp/out.gr", "graph.gr")
    copyfile("tmp/out.info", "out.info")

    log.log("Scaffolds save to " + str(directory) + "/scaffolds.fa")
    log.log("Graph save to " + str(directory) + "/graph.gr")
    log.log("Scaffolds description save to " + str(directory) + "/out.info")
    log.log("Thank you for use IGAR!")
    return

args = parse_args()
run(args)
