# ContigMerger
    Пакет для склейки двух контигов. Есть два контига и выравненые на 
    нах риды. Склеивает два контига в один с каким-то расстоянием 
    между ними. И перевыравнивает риды на новый контиг. Нужен 
    для визуализации выранивания в случе проведения ребра между 
    контигами. 
    
* void evaluate(std::string contigsINFileName, std::string samReads1FileName,
                    std::string samReads2FileName, std::string contigOUTFileName,
                    std::string samOutFileName,
                    std::string contig1Name, std::string contig2Name) 
                    функция, которая склеивает два контига 
                    и перевыранивает риды. 
                     
    * contigsINFileName - файл со всеми контигами в fasta/fastq формате. 
    * samReads1FileName - файл с ридами, выравнеными на первый контиг в SAM/BAM формате. 
    * samReads2FileName - файл с ридами, выравненый на второй контиг в SAM/BAM формате. 
    * contigOUTFileName - файл, в который будет выведен новый скленый контиг в fasta/fastq формате. 
    * samOutFileName - файл, в который будут выведены первый и вторые риды выравненные на новый контиг. 
    * contig1Name - имя первого контига.
    * contig2Name - имя второго контига.
    