## Пакет Builder
### Архитектура
* Изначально управление передается GraphControl, который отвечает за парсинг аргументов и пониманет кому дальше передавать управление. 
* Имеется общий ContigGraph, который и будет строиться по ходу выполнения программы. 
* В зависимости от входных параметров GraphControl передает ContigGraph реализациям абстрактного класса GraphBuilder которые в 
свою очередь достраивают данный граф. Причем каждый следующий GraphBuilder не будет иметь право портить то, что сделал 
предыдущий GraphBuilder. 
* У каждого GraphBuilder имеется SamFileWriteEdge. Он отвечает 
за то, что бы в отдельной папочки для каждого ребра выводилась 
информация о ридах на основание которых данное ребро было проведено. 
* В местах, где используется SystemTools важно, что бы операционая система была Linux и что бы был установлен определенный soft. 

![](../../resources/BuilderClassDiagram.jpg)

### GraphControl
    Точка входа, отвечает за парсинг аргументов и последовательного запуска 
    нужных GraphBuilders с нужными параметрами. 
    
* Имеет одну публичную функцию evaluate(int argc, char **argv). 
Аргументы совпадают с аргументами при запуске программы builder.
* Создает целевой граф, который поочередно изменяют разные GraphBuilder. 
* В конце выводит данный граф в graph.gr
* Работа происходит в папке contig_graph_<текущая дата>

Для добавления очередной библиотеки нужно написать 
-n <тип библиотеки> а далее указать флаги необходимые 
для этой библиотеки. 

Что бы получить полный список возможных флагов можно запусть ./build --help

Пример запуска:
./build -n RNA_SPLIT -l splitLib -r /path/to/contig.fasta -p /path/to/read.fasta
-n RNA_PAIR -l rnaPairLib -f /path/to/first/read.sam -s /path/to/second/read.sam

### Tools
    Пакет, отвечающий за какие-то мелкие вспомогательные команды
#### SeqanUtils
    Несколько функций, которые работают с пакетом Seqan
* static std::string cutReadName(seqan::BamAlignmentRecord read); 
Возвращает имя рида без "/1" и "/2" на конце.  
* static std::string dna5ToString(seqan::Dna5* seq, int len);
Переводит dna5 в строчку
* static void writeRec(seqan::SeqFileOut& out, std::string name, std::string seq);
Выводит запись в fasta/fastq формате
#### SystemAlignmentTools
     Запуск разных видов выравниваний. 
   * Операется, что установлен конкретный soft и работа происходит 
   под UNIX-like ОС. 
   * static void alignmentRNA(std::string refFileName, std::string rnaFileName, std::string resFileName, std::string path=".");
   Выранивает РНК риды на refFile. refFileNAme и rnaFileName файлф в fasta/fastq формате. resFileName файл с выравниванием
   в sam формате. Path путь до папки, где все будет запускаться. Путь без '/' в конце. resFileName уазывать только имя
   файла, а не путь. Файл с таким именем создастся в папке path. refFileNAme и rnaFileName нужно указывать полный или 
   относительный путь. 
   * static void alignmentREF(std::string refFileName, std::string queryFileName) 
   Выравнивает query на ref. queryFileName - имя файла в fasta/fastq формате контигов, которые
   хотим выранить. refFileName - имя файла в fasta/fastq формате с эталлоной сборкой. Результатом
   будет файлик out.coords
   
### SamFileWriteEdge
        Выводит для каждого ребра информацию о ридах из-за которого данное ребро 
        образовалось. 
        
* dir - папка, в которой будут создаваться фалики и в которые будет записываться информация о ридах. 
 В этой папке будет дополнительно создана папка edge и дальнешие записи будут в ней. 
* void setFileIn(seqan::BamFileIn* in) 
необходимо чисто что бы получить конест для создания BamFileOut
* virtual void writeEdge(int edgeID, seqan::BamAlignmentRecord read1, seqan::BamAlignmentRecord read2);
добавляет информацию для ребра edgeID, что оно появилось на основание парных ридов read1 и read2