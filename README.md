# SGTK
1. <a href="#sec1">About SGTK</a><br>
2. <a href="#sec2">Installation</a><br>
3. <a href="#sec3">Running SGTK</a><br>
&nbsp;&nbsp;&nbsp;&nbsp;3.1. <a href="#sec3.1">Input</a><br>
&nbsp;&nbsp;&nbsp;&nbsp;3.2. <a href="#sec3.2">Command line options</a><br>
&nbsp;&nbsp;&nbsp;&nbsp;3.3. <a href="#sec3.3">Output</a><br>
&nbsp;&nbsp;&nbsp;&nbsp;3.4. <a href="#sec3.4">Example</a><br>
4. <a href="#sec4">References</a><br>
5. <a href="#sec5">Feedback and bug reports</a><br>

<a name="sec1"></a>
## 1. About SGTK
SGTK &ndash; Scaffold Graph ToolKit &ndash; is a tool which is capable of scaffold graph construction and
interactive visualization using different kinds of sequencing data. This manual will help you to install and run SGTK.

The current version of SGTK takes as input a set of contigs and an arbitrary number of linkage information sources, such as:
    <li>paired reads</li>
    <li>long reads</li>
    <li>paired and unpaired RNA-seq reads</li>
    <li>scaffolds</li>
    <li>assembly graph in GFA, FASTG formats</li>
    <li>reference sequences</li>

<a name="sec2"></a>
## 2. Installation
To install SGTK download the source code and compile it yourself.

To compile SGTK by yourself you will need the following libraries to be pre-installed:
    <li>python3</li>
    <li>g++ </li>
    <li>cmake </li>
    <li>zlib</li>
    <li>Threads</li>
    <li>Boost</li>
    <li>SEQAN</li>
    <li>minimap2</li>
    <li>STAR</li>
    <li>nucmer</li>


If you meet these requirements, you can download the SGTK source code:

    git clone https://github.com/olga24912/SGTK.git
    cd SGTK/scaffolder/
and build it with the following script:

    ./compile.sh
SGTK will be built in the directory <code>./bin</code>. If you wish to install SGTK into another directory, you can specify the full path of destination folder by running the following command in <code>bash</code> or <code>sh</code>:

    PREFIX=<destination_dir> ./compile.sh
for example:

    PREFIX=/usr/local ./compile.sh

which will install SGTK into <code>/usr/local/bin</code>.

After installation you will get the following files in <code>./bin</code> (or <code>&lt;destination_dir>/bin</code> if you specified PREFIX) directory:
    <li><code>rna_scaffolder.py</code>  (main executable script for building rna scaffolders)</li>
    <li><code>visualize_scaffold_graph.py</code>  (main executable script for visualization scaffold graph) </li>
    <li><code>build</code>  (graph construction module)</li>
    <li><code>filter</code>  (graph simplification and building scaffolds module)</li>
    <li><code>mergeGraph</code>  (graph merging module)</li>
    <li><code>readSplitter</code> (module for splitting RNA-seq reads)</li>
    <li><code>mainPage.html</code> (main html page for visualization)</li>
    <li><code>scripts/</code> (folder containing JS necessary for visualization)</li>
    <li><code>scripts/script.js</code> </li>
    <li><code>scripts/light.css</code> </li>
    <li><code>scripts/search.js</code> </li>
    <li><code>scripts/showListScript.js</code> </li>
    <li><code>scripts/zoomChoose.js</code> </li>
    <li><code>scripts/scaffoldgraph.js</code> </li>
    <li><code>scripts/Icon.png</code> </li>
    <li><code>scripts/defaultWeight.js</code> </li>
    <li><code>scripts/findAreaScript.js</code> </li>
    <li><code>scripts/freeLayout.js</code> </li>
    <li><code>scripts/handleScaffolds.js</code> </li>
    <li><code>scripts/handleAlongChromosomes.js</code> </li>
    <li><code>scripts/handleDiffInLibs.js</code> </li>
    <li><code>scripts/external/cytoscape.js</code></li>
    <li><code>scripts/external/cytoscape-dagre.js</code></li>
    <li><code>scripts/external/dagre.js</code></li>
    <li><code>scripts/external/cytoscape-qtip.js</code></li>
    <li><code>scripts/external/cytoscape.js-navigator/bower.json</code></li>
    <li><code>scripts/external/cytoscape.js-navigator/package.json</code></li>
    <li><code>scripts/external/cytoscape.js-navigator/cytoscape-navigator.js</code></li>
    <li><code>scripts/external/cytoscape.js-navigator/gulpfile.js</code></li>


We also suggest adding SGTK installation directory to the <code>PATH</code> variable.

<a name="sec3"></a>
## 3. Running SGTK
<a name="sec3.1"></a>
### 3.1. Input
SGTK takes as input contigs, reference and an arbitrary number of linkage
information sources.
#### Contigs
Contigs in FASTA format. You can provide a few contigs files, in this case they will be merged together.
Make sure that all contigs have different names.

If you have contigs in forward and reverse-complementary copies, provide only one copy.
#### Reference
If you provide reference genome in FASTA format the genome browser visualization will be available.
In case of providing several files they will be merged together and chromosomes’ names will be changed
depending on the references’ files names.

If you have chromosomes in forward and reverse-complementary copies, provide only one copy.

#### Paired-end and mate-pair reads
SGTK takes paired reads in FASTA or FASTQ format in two files or alignments in SAM/BAM format in two files as well.
In case of providing reads in FASTA/FASTQ format they will be aligned on contigs by minimap2:

    minimap2 -ax sr <contigs_file> <dna1> > dna1.sam
    minimap2 -ax sr <contigs_file> <dna2> > dna2.sam


#### RNA-seq reads
You can provide both paired and single RNA-seq reads. Paired reads will be aligned to the contigs by STAR independently. Single reads will be split into two parts and after that will be aligned by STAR.

#### PacBio/Oxford Nanopore reads
SGTK takes PacBio/Oxford Nanopore reads in FASTA or FASTQ format and aligns them on contigs by minimap2:

    minimap2 -x map-pb <contigs_file> <long_reads> > out.paf


#### Scaffolds
There are two ways to provide scaffolds:
    <li>scaffolds in FASTA format. In this case contigs will be aligned to scaffolds by nucmer.</li>
    <li>scaffolds paths in INFO format </li>

In INFO format each line describes scaffold in the following format:<br>

    >SCAFFOLD_NAME (CONTIG_NAME_0 CONTIG_ID_0 DIRECTION_0) (CONTIG_NAME_1 CONTIG_ID_1 DIRECTION_1)

For example, one line in INFO format:

    >scaffoldName (contigName0 0 +) (contigName1 1 -) (contigName2 2 -) (contigName3 3 +)

#### Assembly graph
Assembly graph can be provide in FASTG or in GFA formats. It will be used instead of contigs and connection from
assembly graph will be detected. In case of graph in GFA format scaffolds(paths) also will be visualized.

#### Connection list
Also it is possible to provide a file with your own connections, where each line describes one connection in the following format:

    (CONTIG_NAME_0 CONTOG_DIRECTION_0) (CONTIG_NAME_1 CONTOG_DIRECTION_1) WEIGHT LEN "EXTRA_INFO"

For example:

    (contig0 -) (contig1 +) 32.5 1168 "some info about connection"

#### Scaffold Graph
All information about connections will be merged in one scaffold graph in internal format. You can provide
extra scaffold graph on your own.

Scaffold graph has the following format:

    LIB_DESCRIPTION
    NODES_DESCRIPTION
    EDGES_DESCRIPTION

The format for LIB_DESCRIPTION:

    LIB_NUMBER
    (l LIB_ID LIB_COLOR LIB_NAME LIB_TYPE)*LIB_NUMBER

Where LIB_TYPE = {CONNECTION | LONG | DNA_PAIR | RNA_PAIR | RNA_SPLIT_50 | RNA_SPLIT_30 | GFA | FASTG | SCAFF}
<br>
<br>
The format for NODES_DESCRIPTION:

    NODES_NUMBER
    (v NODE_ID NODE_NAME NODE_LEN)* NODES_NUMBER

Format for EDGES_DESCRIPTION:

    EDGES_NUMBER
    (e EDGE_ID NODE_ID_0 NODE_ID_1 LIB_ID WEIGHT LEN "EXTRA_INFO")*EDGES_NUMBER

<br>

    1
    l 0 #ff0000 lib_name DNA_PAIR
    4
    v 0 node_0 20132
    v 1 node_0-rev 20132
    v 2 node_1 20400
    v 3 node_1-rev 20400
    2
    e 0 0 2 0 32.5 200 "some extra info"
    e 1 3 1 0 32.5 200

Note, that you should specify contigs, FASTG or GFA file.

<a name="sec3.2"></a>
### 3.2. Command line options
<p>
    To run scaffold graph visualization from the command line, type


    visualize_scaffold_graph.py [options]

Note that we assume that SGTK installation directory is added to the <code>PATH</code> variable (otherwise provide full path to SGTK executable: <code>&lt;installation dir>/visualize_scaffold_graph.py</code>).

#### Options
<p>
    <code>-h</code> (or <code>--help</code>)
    &nbsp;&nbsp;&nbsp;&nbsp;Print help.
</p>
<p>
    <code>-c</code> (or <code>--contig</code>) <code> &lt;file_name> </code>
    &nbsp;&nbsp;&nbsp;&nbsp;File with contigs in fasta format.
</p>
<p>
    <code>-s</code> (or <code>--scaffolds</code>) <code> &lt;file_name> </code>
    &nbsp;&nbsp;&nbsp;&nbsp;File with scaffolds in fasta format.
</p>
<p>
    <code>--fastg &lt;file_name> </code>
        File with assembly graph in FASTG format.
</p>
<p>
    <code>--gfa &lt;file_name> </code>
    File with assembly graph in GFA format.
</p>
<p>
    <code>--fr &lt;file_name_1> &lt;file_name_2> </code>
    &nbsp;&nbsp;&nbsp;&nbsp; File with left reads and file with right reads for paired-end/mate-pair DNA library with forward-reverse orientation.
</p>
<p>
    <code>--rf &lt;file_name_1> &lt;file_name_2> </code>
    &nbsp;&nbsp;&nbsp;&nbsp; File with left reads and file with right reads for paired-end/mate-pair DNA library with reverse-forward orientation.
</p>
<p>
    <code>--ff &lt;file_name_1> &lt;file_name_2> </code>
    &nbsp;&nbsp;&nbsp;&nbsp; File with left reads and file with right reads for paired-end/mate-pair DNA library with forward-forward orientation.
</p>

<p>
    <code>--fr_sam &lt;file_name_1> &lt;file_name_2> </code>
    &nbsp;&nbsp;&nbsp;&nbsp; File with left reads alignment and file with right reads alignment in SAM/BAM format for paired-end/mate-pair DNA library with forward-reverse orientation.
</p>
<p>
    <code>--rf_sam &lt;file_name_1> &lt;file_name_2> </code>
    &nbsp;&nbsp;&nbsp;&nbsp; File with left reads alignment and file with right reads alignment in SAM/BAM format for paired-end/mate-pair DNA library with reverse-forward orientation.
</p>
<p>
    <code>--ff_sam &lt;file_name_1> &lt;file_name_2> </code>
    &nbsp;&nbsp;&nbsp;&nbsp; File with left reads alignment and file with right reads alignment in SAM/BAM format for paired-end/mate-pair DNA library with forward-forward orientation.
</p>
<p>
    <code>--long &lt;file_name> </code>
    &nbsp;&nbsp;&nbsp;&nbsp; File with PacBio/Oxford Nanopore reads.
</p>
<p>
    <code>--scg &lt;file_name> </code>
    &nbsp;&nbsp;&nbsp;&nbsp;File with connection list.
</p>
<p>
    <code>--gr &lt;file_name> </code>
    &nbsp;&nbsp;&nbsp;&nbsp;File with scaffold graph description.
</p>
<p>
    <code>--rna-p &lt;file_name_1> &lt;file_name_2> </code>
    &nbsp;&nbsp;&nbsp;&nbsp; File with left reads and file with right reads for paired-end RNA library.
</p>
<p>
    <code>--rna-s &lt;file_name> </code>
    &nbsp;&nbsp;&nbsp;&nbsp; File for single-read RNA library.
</p>
<p>
    <code>--ref &lt;file_name> </code>
    &nbsp;&nbsp;&nbsp;&nbsp; File with reference in FASTA format.
</p>
<p>
    <code>--scafinfo &lt;file_name> </code>
    &nbsp;&nbsp;&nbsp;&nbsp; File with scaffolds in info format.
</p>
<p>
    <code>--label &lt;label1 label2 ...></code>
    &nbsp;&nbsp;&nbsp;&nbsp; Labels for libraries in given order.
</p>
<p>
    <code>--color &lt;color1 color2 ...></code>
    &nbsp;&nbsp;&nbsp;&nbsp; Colors for libraries in given order.
</p>
<p>
    <code>-o</code> (or <code>--local_output_dir</code>) <code> &lt;output_dir> </code>
    &nbsp;&nbsp;&nbsp;&nbsp;Specify the output directory. The default output directory is  "./"
</p>

<a name="sec3.3"></a>
<h4>3.3. Output</h4>

SGTK stores all output files in <code>&lt;output_dir> </code>, which is set by the user.

<ul>
    <li><code>&lt;output_dir>/scripts/</code> files with scaffold graph description. It is needed for main.html to work properly </li>
    <li><code>&lt;output_dir>/main.html</code> main file for graph visualization </li>
</ul>

Run main.html to see visualisation.

<a name="sec3.4"></a>
### 3.4. Example
To test the toy data set, you can also run the following command from the SGTK <code>bin</code> directory:

    python3 visualize_scaffold_graph.py -c ../share/test_dataset/contigs.fasta --rna-p ../share/test_dataset/read_1.fasta \
    ../share/test_dataset/read_2.fasta -o test

If you have several paired-end reads, scaffolds and reference:
<li> contigs

    ../share/test_dataset/contigs.fasta

<li> paired-end library

    ../share/test_dataset/read_1.fasta
    ../share/test_dataset/read_2.fasta

<li> scaffolds

    scaf.info

<li> reference

    ref.fasta


If you would like to set labels and colors, you need to set labels and colors for all libraries in order of definition

    python3 visualize_scaffold_graph.py -c ../share/test_dataset/contigs.fasta \
    --fr ../share/test_dataset/read_1.fasta ../share/test_dataset/read_2.fasta \
    --scafinfo ../share/test_dataset/scaf.info \
    --ref ../share/test_dataset/ref.fasta \
    --label pair scaf ref \
    --color "#0000ff" "#00ff00" "#ff0000" \
    -o output

<a name="sec4"></a>
<h2>4. Visualization</h2>

![Main screen](https://raw.githubusercontent.com/olga24912/SGTK/develop/resources/pic/main.png)

<a name="sec5"></a>
<h2>5. References </h2>

<a name="sec6"></a>
<h2>6. Feedback and bug reports </h2>