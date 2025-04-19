# Genomica Evolutiva (UNMSM-INS)
Códigos empleados en las lecciones de la Asignatura de "Genómica Evolutiva" del Diplomado en Bioinformática aplicada a la Salud Pública
![INS_UNMSM](https://github.com/user-attachments/assets/f6791aaa-c3ca-4353-93bf-a9b6cc33067c)

## Profesor: Blgo. Victor Jiménez Vásquez
![imagen](https://github.com/user-attachments/assets/89700c0f-916f-401c-b885-b22fc228b7fd)


correos: vr.jimenez.vs@gmail.com, vjimenez@ins.gob.pe, docenteposgrado14.fcb@unmsm.edu.pe

LINKED IN: www.linkedin.com/in/victor-jiménez-vásquez-9bb59666

CTI-VITAE: https://ctivitae.concytec.gob.pe/appDirectorioCTI/VerDatosInvestigador.do?id_investigador=17249

ORCID: https://orcid.org/my-orcid?orcid=0000-0001-6547-6999

GOOGLE: https://scholar.google.com/citations?hl=en

ENTREVISTAS DESTACADAS: https://go.imedia.pe/36QhX, https://go.imedia.pe/3WPZk, https://canaln.pe/actualidad/biologo-ins-sobre-nueva-variante-covid-todavia-no-ha-llegado-al-peru-n464914, https://www.youtube.com/watch?v=HYtZ_JsyHNI&t=7s, https://www.youtube.com/watch?v=M-o-2a2tD8E&t=3s, https://www.youtube.com/watch?v=j1AOTxfqnRw, https://www.youtube.com/watch?v=DjD9Ip8D1_4

# Leccion 1 : Práctica I: Obtención de información de secuenciación genómica.
![Captura desde 2025-04-18 21-27-20](https://github.com/user-attachments/assets/14ddf24b-bac9-4d0b-aaf8-2acc39725656)
```r
###################################
## A. DESCARGA DE ARCHIVOS FASTQ ##
###################################

#paso 1 - ingresar a la página base
https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
#paso 2 - descargar el programa
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz -O stk.tar.gz
#paso 3 - descomprimir el archivo
chmod 777 stk.tar.gz
#paso 4 - enviar el programa a la carpeta bin
export PATH=$PATH:$PWD/sratoolkit.3.0.10-ubuntu64/bin
#paso 5 - descargar archivos de secuenciacion (fastq) a partir de una lista en formato "txt"
prefetch --max-size 50G --option-file sra_accessions_1.txt ;
#paso 6 - explorar los argumentos del comando
prefetch -h
#paso 7 - mover de directorio
mv */*.sra . ;
#paso 8 - eliminar las carpetas iniciales (configurar de acuerdo con el archivo "txt")
rm -r ERR12389866/ ERR12543675/
#paso 9 - extraer los archivos fastq 
fasterq-dump --split-files *.sra ;
#paso 10 - comprimir los archivos fastq
gzip *fastq ;
#paso 11 - estimar la calidad de las lecturas
fastqc *

##########################################################
## B. DESCARGA DE ARCHIVOS ENSAMBLADOS EN FORMATO FASTA ##
##########################################################

#paso 1: instalar NCBI-DATASETS
conda create -n ncbi_datasets
conda activate ncbi_datasets
conda install -c conda-forge ncbi-datasets-cli

#paso 2: emplear el archivo "accessions.txt"

#paso 3: descargar los genomas con la lista sugerida. Emplear el codigo "command_ncbidatasets.sh" disponible en https://github.com/Vjimenez-vasquez/NCBI-DATASETS

./command_ncbidatasets.sh accessions.txt 
```

# Leccion 2 : Práctica II: Ensamblaje y anotación de genomas bacterianos.
![Captura desde 2025-04-18 21-51-53](https://github.com/user-attachments/assets/a4922bb3-c901-4fc8-8342-a4c0a0021681)
```r
#################
## A. TRIMMING ##
#################

# paso 1: descargar TRIMMOMATIC
web oficial: http://www.usadellab.org/cms/?page=trimmomatic
link directo de descarga del "binary": http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip

# paso 2: emplear TRIMMOMATIC (loop)
for r1 in *fastq.gz
do
prefix=$(basename $r1 _L001_R1_001.fastq.gz)
r2=${prefix}_L001_R2_001.fastq.gz
java -jar trimmomatic-0.39.jar PE -threads 28 $r1 $r2 ${prefix}_f_paired.fq.gz ${prefix}_f_unpaired.fq.gz ${prefix}_r_paired.fq.gz ${prefix}_r_unpaired.fq.gz SLIDINGWINDOW:4:20 MINLEN:50 ;
done ; 

#paso 3: mover los FASTQ resultantes a una carpeta 
mkdir trimm ; 
mv *_paired.fq.gz trimm/ ; 
rm *_unpaired.fq.gz ;
cd trimm/ ;

# paso 4: analizar la calidad resultante
fastqc *.gz -t 25 ; 
mkdir fastqc ;
mv *.html *.zip fastqc/ ; 
ls -lh ;

###################
## B. ENSAMBLAJE ##
###################

#paso 1: instalar SPADES
conda create -n spades
conda install bioconda::spades
conda activate spades

#paso 2: ensamblar con SPADES (loop)
for r1 in *fq.gz
do
prefix=$(basename $r1 _f_paired.fq.gz)
r2=${prefix}_r_paired.fq.gz
spades --pe1-1 $r1 --pe1-2 $r2 --careful --cov-cutoff 100 -t 30 -m 15 -o ${prefix}_spades ;
mv ${prefix}_spades/scaffolds.fasta ${prefix}_spades/${prefix}_spades_scaffolds.fasta ;
mv ${prefix}_spades/${prefix}_spades_scaffolds.fasta . ;
done ;
rmdir *fq_spades ;
mv *scaffolds.fasta .. ;

###########################
## C. ANOTACION GENOMICA ##
###########################

# paso 1: instalacion de prokka
conda create -n prokka_env
conda activate prokka_env
conda install -c conda-forge -c biocondaconda install conda-forge::r-base prokka

# paso 2: analisis en prokka (loop)
mkdir annotation/ ;
for r1 in *fa
do
prefix=$(basename $r1 .fa)
prokka --cpus 4 $r1 -o ${prefix} --prefix ${prefix} --kingdom Viruses ; 
mv ${prefix}/*.gff annotation/${prefix}.gff
done ;

# paso 3: instalacion de artemis para visualización
conda create -n art
conda activate art
conda install bioconda::artemis
conda install conda-forge::r-base
```

# Leccion 3 : Práctica I: Identificación de factores de virulencia bacteriana.
![Captura desde 2025-04-19 08-47-42](https://github.com/user-attachments/assets/c5c47bd5-5104-4e67-a821-caabb312068c)
```r
#############
## A. VFDB ##
#############

#paso 1 : http://www.mgc.ac.cn/VFs/
#paso 2 : Default webpage accessible to all users worldwide
#paso 3 : Download
#paso 4 : DNA sequences of full dataset
#paso 5 : Protein sequences of full dataset
#paso 6 : descomprimir
gzip -d VFDB_setB_nt.fas.gz 
gzip -d VFDB_setB_pro.fas.gz

##############
## B. BLAST ##
##############

#paso 1 : instalacion a traves de CONDA
conda install bioconda::blast
or
conda install -c conda-forge -c bioconda -c defaults blast

#paso 2 : correr BLAST
makeblastdb -in VFDB_setB_nt.fas -dbtype nucl ;
blastn -db VFDB_setB_nt.fas -query GCA_001183825.1.fasta -perc_identity 90 -outfmt 6 -num_threads 4 > blast.csv ;
head blast.csv ;
cat blast.csv ;

#paso 3 : añadir headers al resultado
sed '1i query.acc.ver subject.acc.ver perc.identity alignment.length mismatches gap.opens q.start q.end s.start s.end evalue bit.score' blast.csv | tr " " "\t" > blast.2.csv

#paso 4 : revisar resultados
head blast.2.csv
cat blast.2.csv

#######################
## C. Analizar con R ##
#######################

#paso 1 : instalar R
conda install -c conda-forge -c bioconda -c defaults r-base

#paso 2 : leer la data resultante de blast#
data <- read.csv("blast.2.csv", sep="\t", header=TRUE)
#conocer el número de filas y columnas de la tabla resultante#
dim(data)

#paso 3 :conocer las filas asignadas a una columna determinada#
length(data$subject.acc.ver)

#paso 4 : conocer el número de elementos únicos de esa columna#
length(unique(data$subject.acc.ver))
length(unique(data$query.acc.ver))

#paso 5 : conocer estadísticos básicos en un solo paso#
summary(data$query.acc.ver)
summary(data$alignment.length)

#paso 6 : obtener un boxplot de los porcentajes de identidad#
boxplot(data$perc.identity)
boxplot(data$perc.identity, xlab="genoma", ylab="% identidad")
summary(data$perc.identity)
data.frame(names(data))

#paso 7 : obtener un plot longitud de alineamiento vs %identidad#
plot(data$alignment.length, data$perc.identity, xlab="length", ylab="% identity", main="BLASTn VFDB vs Chlamydia", pch=16, col="blue", cex=2)

#paso 8: generar un archivo "bed" con R
seq <- data.frame(genome=data$query.acc.ver, start=data$q.start, end=data$q.end)
head(seq)
write.table(seq, "extract.txt", sep="\t", row.names = F, col.names =F, quot=F)

#################
## D. BEDTOOLS ##
#################

#paso 1 : instalar bedtools desde conda para extraer las regiones "blasteadas"
conda install conda install -c conda-forge -c bioconda -c defaults bedtools

#paso 2 : extraer los VFs en formato FASTA
bedtools getfasta -fi  GCA_001183825.1.fasta -bed extract.txt -fo virulence.fasta

#paso 3 : extraer información de los headers de VFDB 
grep ">" VFDB_setB_nt.fas | sed -e 's/]\ \[/*/g' | sed -e 's/]//g' | sed -e 's/\ \[/*/g'| sed -e 's/)\ /*/g' | sed -e 's/*(/*/g' | head -n 10 > headers.txt

#paso 4 : traducir las secuencias con virtual ribosome
https://services.healthtech.dtu.dk/services/VirtualRibosome-2.0/
```

# Leccion 4 : Práctica II: Ensamblaje de genomas bacterianos a partir de secuencias nanopore y evaluación de la calidad.
![Captura desde 2025-04-19 08-50-35](https://github.com/user-attachments/assets/7d9bc79d-bf1c-44e7-b56f-df8bfea3a77f)
```r
## material de apoyo > https://denbi-nanopore-training-course.readthedocs.io/en/stable/index.html ##

#################################
## A. INSTALACION DE PROGRAMAS ##
#################################
#paso 1 : NanoPlot - calidad de secuencias Nanopore
conda install -c conda-forge -c bioconda nanoplot

or

pip install NanoPlot
pip install NanoPlot --upgrade

#paso 2 : Nanofilt : Filtrado por calidad de lecturas Nanopore
conda install -c conda-forge -c defaults -c bioconda nanofilt

or 

pip install nanofilt
pip install nanofilt --upgrade

#paso 3 : Flye: de-novo assembly
conda install -c bioconda flye

or 

git clone https://github.com/fenderglass/Flye
cd Flye
python setup.py install

#paso 4 : Minimap2 - polishing (parte 1)
conda install -c bioconda minimap2

or

git clone https://github.com/lh3/minimap2
cd minimap2 && make

#paso 5 : Racon - polishing (parte 2)
conda install -c bioconda racon

or 

git clone --recursive https://github.com/lbcb-sci/racon.git racon
cd racon
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd build/bin/ 
export PATH=$PATH:$HOME/bin
cp racon $HOME/bin
chmod +x $HOME/bin/racon

#paso 6 : Requerimientos de MEDAKA (Pyabpoa, bcftools, samtools (v1.11), minimap2)
pip install pyabpoa
sudo apt install bcftools
conda install -c bioconda samtools==1.11

#paso 7 : MEDAKA, secuencias consenso (si MEDAKA no funciona correctamente, instala los programas requeridos)
conda install -c conda-forge –c bioconda medaka

or

pip install medaka

###################
## B. ENSAMBLAJE ##
###################

#paso 1 : descargar la informacion (códigos SRR17110067 y SRR17110070)
mkdir sra_files ;
prefetch --max-size 50G --option-file accessions.txt ;
mv */*.sra . ;
fasterq-dump --split-files *.sra 
gzip *.fastq ;
mkdir sra_files ;
mv *.sra sra_files/ ;

#paso 2 : analizar longitudes
zcat SRR17110067.fastq.gz | grep -n "length" | cut -f2 -d'=' | sort -r -n | uniq | head -n 20
zcat SRR17110070.fastq.gz | grep -n "length" | cut -f2 -d'=' | sort -r -n | uniq | head -n 20

#paso 3 : NanoPlot
NanoPlot -t 2 -o SRR17110067_QC --fastq SRR17110067.fastq.gz
NanoPlot -t 2 -o SRR17110070_QC --fastq SRR17110070.fastq.gz

#paso 4 : NanoFilt
gunzip -c SRR17110067.fastq.gz | NanoFilt --logfile nanofilt.log -l 500 -q 10 | gzip > SRR17110067.trim.fastq.gz ;
gunzip -c SRR17110070.fastq.gz | NanoFilt --logfile nanofilt.log -l 500 -q 10 | gzip > SRR17110070.trim.fastq.gz ;
ls -lh ;

#paso 5 : Flye
flye -o SRR17110067.genoma --nano-raw SRR17110067.trim.fastq.gz --threads 4 ;
flye -o SRR17110070.genoma --nano-raw SRR17110070.trim.fastq.gz --threads 4 ;
ls -lh ;

#paso 6 : Minimap2 + Racon (Polishing)
minimap2 -x ava-ont -t 4 SRR17110067.genoma/assembly.fasta SRR17110067.trim.fastq.gz > overlaps1.paf ;
racon -t 4 SRR17110067.trim.fastq.gz overlaps1.paf SRR17110067.genoma/assembly.fasta > SRR17110067.racon1.fasta ;

minimap2 -x ava-ont -t 4 SRR17110070.genoma/assembly.fasta SRR17110070.trim.fastq.gz > overlaps2.paf ;
racon -t 4 SRR17110070.trim.fastq.gz overlaps2.paf SRR17110070.genoma/assembly.fasta > SRR17110070.racon1.fasta ;

minimap2 -x ava-ont -t 4 SRR17110067.racon1.fasta SRR17110067.trim.fastq.gz > overlaps3.paf ;
racon -t 4 SRR17110067.trim.fastq.gz overlaps3.paf SRR17110067.racon1.fasta > SRR17110067.racon2.fasta ;

minimap2 -x ava-ont -t 4 SRR17110070.racon1.fasta SRR17110070.trim.fastq.gz > overlaps4.paf ;
racon -t 4 SRR17110070.trim.fastq.gz overlaps4.paf SRR17110070.racon1.fasta > SRR17110070.racon2.fasta ;

#paso 7 : Medaka (consensus)
medaka_consensus -i SRR17110070.trim.fastq.gz -d SRR17110070.racon2.fasta -o medaka_SRR17110070 -t 4 ;
medaka_consensus -i SRR17110067.trim.fastq.gz -d SRR17110067.racon2.fasta -o medaka_SRR17110067 -t 4 ;

#paso 8 : QUAST
quast.py -o quast_results -m 0 consensus.fasta

#paso 9 :Bandage
```

# Leccion 5 : Práctica I: estrategias para la obtención de árboles de máxima verosimilitud.
![Captura desde 2025-04-19 09-59-39](https://github.com/user-attachments/assets/6f9dc036-ac38-457a-93c6-1df5b65691d1)
```r

#paso 1 : alineamiento
mafft --add output.fasta --nomemsave --keeplength --thread 4 --retree 1 --adjustdirection --reorder reference.fasta > output3.fasta ;
aliview output3.fasta ;
ls -lh

#paso 2 : comandos para obtener un arbol ML
raxmlHPC-PTHREADS -p 123568 -m GTRCAT -s sequences.phy -T 30 -# 4 -n nwk ; 
rm raw_tree.nwk ;
mv RAxML_bestTree.nwk raw_tree.nwk ;
mkdir raxml ; 
mv RAx* raxml/ ;

```

# Leccion 6 : Práctica II: Compilación de metadatos.
![Captura desde 2025-04-19 11-19-42](https://github.com/user-attachments/assets/93b0ccf6-bf54-40c1-b1f2-60e974875dc9)
```r
```

# Leccion 7 : Práctica I: Modelos demográficos.
![Captura desde 2025-04-19 11-26-39](https://github.com/user-attachments/assets/587c3511-d121-4e7c-a90f-049943530870)
```r
####################################################
## A. INSTALAR APLICACIONES Y PROGRAMAS (WINDOWS) ##
####################################################
#paso 1: web oficial de BEAST
https://beast.community/
https://beast.community/install_on_windows
https://github.com/beast-dev/beast-mcmc/releases
BEAST v1.10.4
BEAST.v1.10.4.zip

#paso 2: web oficial de JAVA
Installing JAVA
https://www.oracle.com/java/technologies/downloads/?er=221886
Windows
x64 Installer
https://download.oracle.com/java/24/latest/jdk-24_windows-x64_bin.exe (sha256)

#paso 3 : descargar BEAGLE
https://beast.community/beagle
https://github.com/beagle-dev/beagle-lib
BEAGLE v4.0.0 for Windows 64-bit

#paso 4 : descargar TRACER
https://github.com/beast-dev/tracer/releases/tag/v1.7.2

##################################################
## B. genera el archivo "XML" en BEAUti v1.10.4 ##
##################################################

#paso 1 :
File
Import Data...
Cargar el alineamiento FASTA, observar el resultado en la pestaña "Partitions"

#paso 2 :
Click en la pestaña "Tips"
check en "Use tip dates"
click en "Import Dates"
seleccionar "Parse as a number" o "Parse as a calendar date" según sea el caso
Tip date sampling (Sampling uniformly from precision)

#paso 3 :
click en la pestaña "Sites"
Substitution Model "HKY"

#paso 4 :
click en la pestaña "Clocks"
clock Type "Uncorrelated relaxed clock"

#paso 5 :
click en la pestaña "Trees"
Tree Prior: "Coalescent: Bayesian Skyline"

#paso 6 :
click en la pestaña "Priors"
click + ENTER en todas las opciones

#paso 7 :
click en la pestaña "MCMC"
Length of chain: 10000000
Echo state to screen every: 1000
Log parameters every : 1000
click en "Add .txt suffix"
click "Create tree log file with branch length in substitutions"
click en "Generate BEAST File..."
nominar y guardar el archivo "xml"

##################################################################
## C. generar las cadenas de MARKOV-MONTECARLO en BEAST v1.10.4 ##
##################################################################

#paso 1 :
Choose file ...
agregar el archivo xml
run

######################################
## D. Observar resultados en TRACER ##
######################################

#paso 1:
File
Import Trace File
Evaluar que los valores ESS sean mayores a 100 e idealmente mayores a 200
click en "age (root)"
Revisar las pestañas "Estimates" y "Marginal Density"
Seleccionar el TraceFile
click en Analysis
click en Bayesian Skyline Reconstruction...
Trees log File: cargar el archivo con terminacion ".(time).trees"
Age of youngest tip: indicar la fecha de la muestra más antigua

```

# Leccion 8 : Práctica II: Filodinámica y obtención de árboles filogenómicos anotados.
![Captura desde 2025-04-19 09-49-55](https://github.com/user-attachments/assets/62893de4-b2b7-441c-a4e3-88820a15d8bd)

```r
# paso 1 : instalacion

installation : from https://docs.nextstrain.org/en/latest/install.html

conda create -n nextstrain \
      --override-channels --strict-channel-priority \
      -c conda-forge -c bioconda --yes \
      augur auspice nextclade \
      snakemake git epiweeks pangolin pangolearn \
      ncbi-datasets-cli csvtk seqkit tsv-utils
      
conda activate nextstrain ;
nextstrain setup --set-default ambient ; 

#paso 2 : alineamiento

mafft --add output.fasta --nomemsave --keeplength --thread 4 --retree 1 --adjustdirection --reorder reference.fasta > output3.fasta ;
aliview output3.fasta ;
ls -lh

#paso 3 : augur - nextstrain
augur tree -a input.fasta --method iqtree --output raw_tree.nwk --substitution-model GTR --nthreads 31 ; 
augur refine --alignment input.fasta --tree raw_tree.nwk --metadata metadata.tsv --output-tree refine_tree.nwk --output-node-data node_Data.json --timetree --coalescent skyline --gen-per-year 52 --root best --covariance --date-confidence --date-inference marginal --branch-length-inference marginal --year-bounds 1998 2024 --divergence-units mutations-per-site --seed 12548746 ; 
augur ancestral --tree refine_tree.nwk --alignment input.fasta --output-node-data ancestral.json --inference marginal --keep-overhangs ; 
augur translate --tree refine_tree.nwk --ancestral-sequences ancestral.json --reference-sequence sequence.2.gb --output-node-data translate.json ; 
augur traits --tree refine_tree.nwk --metadata metadata.tsv --columns country province lineage result informe isolate year --confidence --output-node-data traits.json ; 
augur export v2 --tree refine_tree.nwk --node-data ancestral.json node_Data.json traits.json translate.json --output denv3.json --auspice-config auspice_config.json --geo-resolutions province country --color-by-metadata country province lineage result informe isolate year --panels tree map entropy frequencies --metadata metadata.tsv --lat-longs lat_longs.tsv ; 
ls -lh ; 

```

# Leccion 9 : Práctica I: Identificación metagenómica y visualización de resultados.
![sankey](https://github.com/user-attachments/assets/fa89c820-36d4-4b4b-9306-274bc7c8d347)

Kraken2 databases availabe at : https://benlangmead.github.io/aws-indexes/k2

Bracken is availabe at : https://github.com/jenniferlu717/Bracken

```r
# path to BRACKEN: /media/ins-bio/DATA01/data_base_download/Bracken-2.7/./bracken
# path to KRAKEN VIRUS DATABASE: /media/ins-bio/DATA01/data_base_download/KRAKENVIRDB
# path to KRAKEN STANDARD DATABASE:/media/ins-bio/DATA01/data_base_download/kraken_standard

# path to KRAKEN BIG-DATABASE (104) : /home/administrador/Documentos/KRAKENPlusDB
# path to KRAKEN16Gb (104) : /home/administrador/Documentos/KRAKEN16Gb
# path to KRAKEN VIRUS DB (104) :  /home/administrador/Documentos/KRAKENVIRDB
# path to BRACKEN (104) : /home/administrador/Documentos/Bracken-2.7

#paso 1 : fastqc
fastqc -t 25 *
mkdir fastqc ; 
mv *.html *.zip fastqc/ ; 
ls -lh ; 

#paso 2 : kraken viral
for r1 in *fastq.gz
do
prefix=$(basename $r1 _L001_R1_001.fastq.gz)
r2=${prefix}_L001_R2_001.fastq.gz
kraken2 --paired --use-names --gzip-compressed --db /media/ins-bio/DATA01/data_base_download/KRAKENVIRDB/ --threads 28 $r1 $r2 --report ${prefix}_report.txt --output ${prefix}_kraken2.out ;
done ;
rm *.fastq.gz_report.txt ; 
mkdir kraken_out ;
mv *.out  kraken_out/ ;
mkdir kraken_txt ; 
mv *.txt kraken_txt/ ;  
cd kraken_txt/ ; 
ls -lh ; 

#paso 3 : bracken
for r1 in *_report.txt
do
prefix=$(basename $r1 _report.txt)
/media/ins-bio/DATA01/data_base_download/Bracken-2.7/./bracken -d /media/ins-bio/DATA01/data_base_download/KRAKENVIRDB/ -i $r1 -o ${prefix}.bracken.txt -l S ; 
done ; 
mkdir species_report ; 
mv *_species.txt species_report/ ;
mkdir bracken_species_abundances ;  
mv /media/ins-bio/DATA01/data_base_download/Bracken-2.7/*.txt bracken_species_abundances/ ; 
mv *.bracken.txt bracken_species_abundances/ ;
ls ;

#paso 4 : PAVIAN en R
if (!require(remotes)) { install.packages("remotes") }
remotes::install_github("fbreitwieser/pavian")
pavian::runApp(port=5000)
```

# Leccion 10 : Práctica II: Identificación de regiones homólogas y análisis pangenómico.
![Captura desde 2025-04-19 09-07-03](https://github.com/user-attachments/assets/a0db9b78-01e4-4024-83fa-42441a60020f)
```r
#################################
## A. INSTALACION DE PROGRAMAS ##
#################################

#paso 1 : instalacion 1
conda install -c conda-forge -c defaults -c bioconda roary
conda install -c conda-forge -c defaults -c bioconda snp-sites
conda install -c conda-forge -c defaults -c bioconda raxml
conda install -c conda-forge -c defaults -c bioconda figtree

#paso 2 : fix roary install problems with sudo
sudo add-apt-repository ppa:c2d4u.team/c2d4u4.0+
sudo apt update
sudo apt-get install roary

###################
## B. ANNOTATION ##
###################
#paso 1 : annotation (PROKKA)

conda activate prokka_env
mkdir -p annotation ;
mkdir -p ffn ;
for r1 in *fasta
do
prefix=$(basename $r1 .fasta)
prokka --cpus 4 $r1 -o ${prefix} --prefix ${prefix} --kingdom Bacteria ; 
mv ${prefix}/*.gff annotation/${prefix}.gff
done ;
conda deactivate ;
cp */*.ffn ffn/ ; 
ls ;

#############################################
## C. IDENTIFICACION DE REGIONES HOMOLOGAS ##
#############################################

#paso 1 : inferring clusters, core genes and accesory genes (ROARY)
# https://github.com/sanger-pathogens/Roary #
roary -p 4 -f roary_output -g 200000 -z -r -e -n -v -cd 80 -i 90 annotation/*.gff
roary -p 4 -f roary_output -g 200000 -r -e -n -v -cd 80 -i 90 annotation/*.gff
cp roary_output/core_gene_alignment.aln . ;
ls -lh ;

##################
## D. FILOGENIA ##
##################
#paso 1 : SNPs alignment (SNP-SITES)
snp-sites -m -o snp1.phy core_gene_alignment.aln ; 
snp-sites -m -c -o snp2.phy core_gene_alignment.aln ; 
ls -lh ;

#paso 2 : phylogeny (RAXML)
raxmlHPC-PTHREADS -p 1254512 -m GTRCAT -s snp2.phy -n nwk -# 20 -T 4 ;
mv RAxML_bestTree.nwk raw_tree.nwk ;
rm RAxML_* ;
mkdir phylogeny ;
mv snp1.phy snp2.phy snp2.phy.reduced raw_tree.nwk core_gene_alignment.aln phylogeny/ ;

#######################################################################
## E. GENERAR UNA TABLA INFORMATIVA DE DATOS ADJUNTOS, METADATA EN R ##
#######################################################################

#paso 1 : cargar el programa "pangenome_command_2.R" en R o R-Studio

#paso 2 : ingresar la ruta correcta en cada caso (donde se encuentran los archivos "gene_presence_absence.csv" y "metadata_1.tsv")
setwd("")
dir()

#paso 3 : ingresar la ruta correcta hasta donde se encuentra el archivo pangenome_command_2.R
source("../../pangenome_command_2.R")

pangenome <- pres_abs(metadata = "metadata_1.tsv", roary_output = "gene_presence_absence.csv", last_column = "3", output = "out_5.tsv")
head(pangenome)
class(pangenome)
pangenome[,1:10]

###################
## F. MICROREACT ##
###################
#paso 1 : visualizacion en microreact
https://microreact.org/
cargar el arbol enraizado (formato .nwk) y la metadata final (out_5.tsv)

```
