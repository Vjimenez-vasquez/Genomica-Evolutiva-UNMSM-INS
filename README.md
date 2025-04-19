# Genomica Evolutiva (UNMSM-INS)
Codigos empleados en las lecciones de la Asignatura de "Genomica Evolutiva" del Diplomado en Bioinformática aplicada a la Salud Pública

## Profesor: Blgo. Victor Jiménez Vásquez

correos: vr.jimenez.vs@gmail.com, vjimenez@ins.gob.pe, docenteposgrado14.fcb@unmsm.edu.pe

CTI-VITAE: https://ctivitae.concytec.gob.pe/appDirectorioCTI/VerDatosInvestigador.do?id_investigador=17249

ORCID: https://orcid.org/my-orcid?orcid=0000-0001-6547-6999

GOOGLE: https://scholar.google.com/citations?hl=en

ENTREVISTAS DESTACADAS: https://go.imedia.pe/36QhX, https://go.imedia.pe/3WPZk, https://canaln.pe/actualidad/biologo-ins-sobre-nueva-variante-covid-todavia-no-ha-llegado-al-peru-n464914, https://www.youtube.com/watch?v=HYtZ_JsyHNI&t=7s, https://www.youtube.com/watch?v=M-o-2a2tD8E&t=3s, https://www.youtube.com/watch?v=j1AOTxfqnRw, https://www.youtube.com/watch?v=DjD9Ip8D1_4

# Leccion 1 : Práctica I: Obtención de información de secuenciación genómica.
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
```r
```

# Leccion 4 : Práctica II: Ensamblaje de genomas bacterianos a partir de secuencias nanopore y evaluación de la calidad.
```r
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
