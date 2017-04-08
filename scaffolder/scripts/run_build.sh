#! /bin/bash
#! /bin/sh

softDir=.
contig=./contigs.fasta
reads1=./reads_1.fastq
reads2=./reads_2.fastq
threadN=1
outputDir=.

while getopts ":c:d:f:s:t:o:" opt; do
    case $opt in
        c)
            contig=$OPTARG
            ;;
        d)
            softDir=$OPTARG
            ;;
        f)
            reads1=$OPTARG
            ;;
        s)
            reads2=$OPTARG
            ;;
        t)
            threadN=$OPTARG
            ;;
        o)
            outputDir=$OPTARG
            ;;
    esac
done

cd $outputDir

lineCnt=(`wc $reads1`)
$softDir/splitFast $reads1 $outputDir/reads1_ $threadN $lineCnt
echo finished split first file
lineCnt=(`wc $reads2`)
$softDir/splitFast $reads2 $outputDir/reads2_ $threadN $lineCnt
echo finished split second file

mkdir genomeDir
$softDir/STAR --runMode genomeGenerate --limitGenomeGenerateRAM 90000000000 --genomeDir genomeDir --runThreadN 20 --genomeFastaFiles $contig

for i in $(seq 0 `expr $threadN - 1`); do
    mkdir gr$i
    cd gr$i
    mv ../reads1_$i ./reads1_$i
    mv ../reads2_$i ./reads2_$i

    $softDir/STAR --runThreadN $threadN --genomeDir ../genomeDir --readFilesIn $reads1
    mv Aligned.out.sam RNA_1.sam
    $softDir/STAR --runThreadN $threadN --genomeDir ../genomeDir --readFilesIn $reads2
    mv Aligned.out.sam RNA_2.sam

    cd $outputDir;
done

for i in $(seq 0 `expr $threadN - 1`); do
    cd gr$i
    cp $softDir/STAR ./STAR
    $softDir/build -n RNA_SPLIT -r $contig -p $reads1 -l split1 -n RNA_SPLIT -r $contig -p $reads2 -l split2 -n RNA_PAIR -f RNA_1.sam -s RNA_2.sam -l pair&
    cd $outputDir;
done

wait

echo -e '0\n0\n0\n' > graphAllLib.gr
for i in $(seq 0 `expr $threadN - 1`); do
    $softDir/mergeGraph ./gr$i/graph.gr ./graphAllLib.gr ./graphAllLib.gr
done

echo 'uploadGraph graph.gr' > filter_config
for i in $(seq 1 `expr $threadN - 1`); do
    echo "mergeLib 0  $(($i * 5)) split1-50-50" >> filter_config
    echo "mergeLib 1  $(($i * 5 + 1)) split1-long-short" >> filter_config
    echo "mergeLib 2  $(($i * 5 + 2)) split2-50-50" >> filter_config
    echo "mergeLib 3  $(($i * 5 + 3)) split2-long-short" >> filter_config
    echo "mergeLib 4  $(($i * 5 + 4)) pair" >> filter_config
done
echo 'print graph.gr' >> filter_config
echo 'exit' >> filter_config

$softDir/filter
