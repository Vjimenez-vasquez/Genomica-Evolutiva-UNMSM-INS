# Genomica Evolutiva (UNMSM-INS)
Codigos empleados en las lecciones de la Asignatura de "Genomica Evolutiva" del Diplomado en Bioinformática aplicada a la Salud Pública

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
## A. DESCARGA DE ARCHIVOS FASTQ

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

## B. DESCARGA DE ARCHIVOS ENSAMBLADOS EN FORMATO FASTA

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
## A. TRIMMING

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

## B. ENSAMBLAJE

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

## C. ANOTACION GENOMICA

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
## A. VFDB

#paso 1 : http://www.mgc.ac.cn/VFs/
#paso 2 : Default webpage accessible to all users worldwide
#paso 3 : Download
#paso 4 : DNA sequences of full dataset
#paso 5 : Protein sequences of full dataset
#paso 6 : descomprimir
gzip -d VFDB_setB_nt.fas.gz 
gzip -d VFDB_setB_pro.fas.gz

## B. BLAST

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

## C. Analizar con R

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

## D. BEDTOOLS

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

## A. INSTALACION DE PROGRAMAS

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

## B. ENSAMBLAJE

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
```r
```

# Leccion 6 : Práctica II: Compilación de metadatos.
```r
```

# Leccion 7 : Práctica I: Modelos demográficos.
```r
```

# Leccion 8 : Práctica II: Filodinámica y obtención de árboles filogenómicos anotados.
```r
```

# Leccion 9 : Práctica I: Identificación metagenómica y visualización de resultados.
```r
```

# Leccion 10 : Práctica II: Identificación de regiones homólogas y análisis pangenómico.
```r
```
