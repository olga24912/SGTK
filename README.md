# bio_scaffolder
Инструмент для построения скаффолдов по ридам РНК и визуализации возможных связей между контигами. 

## Установка

## Запуск
Есть несколько исполняемых файлов. 
  * build --- строит граф связей между контигами. 
  * filter --- фильтрует и визуализирует граф. 
  * addInfoToGraph --- возможность добавить в граф информацию о том, что было соединено в скаффолды. 
  * mergeGraph --- склеить два графа в один. 

### Пример запуска смотрите в папке example

### build
Построение графа связей. 

Вероятный вариант запуска:
    
    ./build -n RNA_PAIR -f RNA1.sam -s RNA2.sam -l rna_pair -n DNA_PAIR -f DNA1.sam -s DNA2.sam -l dna_pair

Каждая отдельная библиотека, на основание которой будут строиться связи начинается 
с флага -n <вид библиотеки>, далее для каждый библиотеки идут свои парметры. 

#### параметры
 * -n <RNA_PAIR|RNA_SPLIT|DNA_PAIR|REF> --- начало новой библиотеке, к данной библиотеке будут относиться все параметры от данного места и до следующего -n или до конца строки.
##### общие параметры на все библиотеки
 * -o --- выводить дополнительную информацию о выравниваниях, на основание которых были проведены связи. 
 Для каждой библиотеки будет создана папка с ее названием и в отдельные файлы с номером ребра будут 
 выписаны все выравнивания на основание которых данное ребро было проведено. 
##### -n RNA_PAIR
Если строим связи на основание парных ридов РНК. 
 *  -f <file with alignment in SAM format>  --- файл с выравниванием первых ридов РНК из пары в sam формате. 
 *  -s <file with alignment in SAM format>  --- файл с выравнивание вторых ридов РНК из пары в sam формате. 
 *  -l <lib name> --- имя библиотеки. 

Как именно можно получить необходимые выравнивания смотрите в примерах. 

##### -n DNA_PAIR
Построение связей по парным ридам ДНК.
 *  -f <file with alignment in SAM format>  --- файл с выравниванием первых ридов РНК из пары в sam формате. 
 *  -s <file with alignment in SAM format>  --- файл с выравнивание вторых ридов РНК из пары в sam формате. 
 *  -l <lib name> --- имя библиотеки. 

Необязательные параметры:
 *  -d <maximum dist between pair reads> ---  максимальная возможная длина между парными ридами. Если длина между парными ридами заведомо больше указанной, то такая пара будет игнорироваться. 
 
##### -n RNA_SPLIT
Построение связей на основание ридов, которые попали на стыки экзонов. 

Для использования этой опции необходимо наличие STAR в папке из 
которой происходит запуск. 

 * -r <contigs file name> --- путь до файла с контигоми между которыми будет осуществляться построение связей. В fasta или fastq формате. 
 * -p <reads file name> --- путь до файла с ридами РНК в fasta или fastq формате. 
 * -l <lib name> --- имя библиотеки. 
 
##### -n REF
Построение связей между контигами на основание уже готовой сборки. 
Есть два варианта для использования этой опции, можно выбрать любой из них 

###### первый вариант
В данном случае необходимо, что бы был установлен mummer (используется nucmer и show-coords)
 * -r <ref file name> --- путь до файла с эталонной сборкой. 
 * -q <contigs file name> --- путь до файла с контигами, которые будут выравниваться на референс
 * -l <lib name> --- имя библиотеки. 
 
Необязательные параметры:
 * -с <min contig len> --- минимальная длина контигов, которые будут учитываться, по умолчанию 500. 

###### второй вариант
 * -t <tsv file name> --- путь до tsv файла с выравниванием. 
 * -l <lib name> --- имя библиотеки

Необязательные параметры:
 * -с <min contig len> --- минимальная длина контигов, которые будут учитываться, по умолчанию 500. 

#### выходные файлы
##### graph.gr
Основной интересующий нас файл с графом. Этот файл должен появиться 
в той же папке, из которой запускалась программа. 

###### структура файла
<количество библиотек - ln> далее идут ln строчек следующего типа: 

l <номер библиотеки> <цвет> <имя библиотеки>  

<количество вершин - vn> далее идут vn строчек следующего типа:  

v <номер вершины> <название контига> <длина контига>  

<количество ребер - en> далее идут en строчек следующего типа:  

e <номер ребра> <вершина, из которой исходит ребро> <вершина, в которую входит ребро> <библиатека, которой принадлежит 
ребро> <вес ребра>  
  
Требуется, что бы номер ребра/библиотеки совпадал с номером строчки в блоке  

###### Пример:  
2  
l 0 #67c669 lib1  
l 1 #ff0000 lib2  
3  
v 0 contig1 1051  
v 1 contig2 150  
v 2 contig3 3800  
2  
e 0 1 2 0 10  
e 1 2 1 1 15
