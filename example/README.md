# Пример
## Файлы
Небольшие тестовые данные, в которые входят 
  * contigs.fasta --- файл с контигами. 
  * ref.fasta  ---  полностью собранный файл.
  * read_1.fasta --- первые риды из пары
  * read_2.fasta --- вторые риды из пары. 

## Подготовка данных
Для начала нам нужно получить файлы RNA1.sam и RNA2.sam с выравнеными 
ридами на контиги. 

Что бы их получить можно, например, использовать STAR. 
    
    mkdir genomeDir
    STAR --runMode genomeGenerate --genomeDir genomeDir --runThreadN 20 --genomeFastaFiles contigs.fasta
    STAR --runThreadN 20 --genomeDir genomeDir --readFilesIn read_1.fasta
    mv Aligned.out.sam RNA1.sam
    STAR --runThreadN 20 --genomeDir genomeDir --readFilesIn read_2.fasta
    mv Aligned.out.sam RNA2.sam
    
    
## Запуск build

   
    ./build -n RNA_PAIR -f RNA1.sam -s RNA2.sam -l rna_pair 
    -n RNA_SPLIT -p read_1.fasta -r contigs.fasta -l rna_split
    -n REF -r ref.fasta -q contigs.fasta -l ref
     
     
На выходе должен получиться файл graph.gr

## Построение scaffold

      
