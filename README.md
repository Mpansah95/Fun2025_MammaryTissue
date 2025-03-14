# MAMMARY TISSUE - BODY SIZE
## FUNCTIONAL GENOMIC (BIOL 6850/5850)

# Group Members
1. Mikailie Caulder
2. Aubrey Wright
3. Lexie Leggett
4. Melika Ghasemi Shiran
5. Prince Mensah Ansah
6. Thais Ribeiro
7. Yagya Adhikari
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Title
The contributions of the Insulin and Insulin-like Signaling (IIS) molecular network in dog body size evolution: contrasting gene expression and function of protein-coding changes

# Background

# Hypotheses
1. The gene expression of the growth hormone/insulin-like signaling molecular network will be significantly different between large and small dog breeds
2. The growth hormone/insulin-like signaling network will contain functional protein-coding sequence variation in the IIS Molecular Pathway that correlates with body size among dogs

# Bioinformatics Plan

## Source of Data
This project will analyze publicly available RNAseq data in the [NCBI SRA database](https://www.ncbi.nlm.nih.gov/sra/?term=). The data that qualifies for our inclusion criteria will be downloaded onto the [Alabama Supercomputer](https://hpcdocs.asc.edu/hpc-tools-categories/sequence-alignment for bioinformatic analysis. Because we are considering differential gene expression analysis in germline tissues - RNAseq data that we sequenced from [mammary gland](https://en.wikipedia.org/wiki/mammary-gland) will be used. Not only that, the organism must be a 
[non-avian reptile](https://www.reddit.com/r/askscience/comments/b5h55n/whats_the_point_of_the_term_nonavian_reptiles). The following data listed below met our criteria stated above. The project ID in [NCBI](https://www.ncbi.nlm.nih.gov/) [BioProject database](https://www.ncbi.nlm.nih.gov/bioproject/) has been listed below


## Data Download

### Downloading fastq files
To test our hypothesis we will use gonad transcriptome data from the following species: Mauremys reevesii, Pogona vitticeps (Georges et al. 2015), Mauremys mutica (Yuan et al. 2021), Anolis carolinensis (Zhang et al. 2017), Trachemys scripta elegans (Hatkevich et al. 2025), and Pelodiscus sinensis (Zhu et al. 2022). Data will be downloaded from the SRA database with the command line tool [SRAtoolkits (v.3.0.0)](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) from NCBI. The fasterq-dump command will be used with the –split-files flag. An accession list, which is a text file containing all the SRA accession numbers for a particular bioproject ID, will be the input file, and the output file will be .fastq files downloaded from the database. To ensure that all the data were successfully downloaded, the file size of the downloaded .fastq file will be manually compared to the file sizes on the NCBI database. For paired-end sequences, the file sizes for R1 and R2 files are expected to be the same. This is a quality step to ensure the reads were completely downloaded.

### Downloading reference genomes and annotation files
Reference genomes for all the organisms included will be downloaded from NCBI through ftp file transfer with the wget Linux command. Gene annotation files (.gff) corresponding to these genomes will also be downloaded. Files will be renamed to ensure easy handling of genomes and corresponding annotation files. Also, each genome and its annotation file will be downloaded into separate directories to avoid downstream errors when processing them.

### Quality Control
[FastQC (v.0.12.1)](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (Andrews, S. 2010) on the Alabama supercomputer will be used to assess the sequence quality control. [MultiQC (v.1.15)](https://github.com/MultiQC/MultiQC) (Ewels et al. 2016) will also be run to summarize all the QC parameters into one file for easy comparison between different sample groups. Sequences that do not pass the quality control test will be trimmed, using [Trimmomatic (v.0.39)](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) (Bolger et al. 2014), to improve base call quality, remove adapters, overrepresented sequences, etc. To ensure that all the trimmed sequences are of good quality, FastQC and MultiQC will be rerun on the trimmed sequences to assess the outcome of the trimming. If all sequences pass the quality check, the trimmed .fastq files will proceed to mapping. Otherwise, another trimming cycle will be done taking into consideration which area of the sequence reads need extra cleaning.

## Mapping

### Genome indexing
[Hisat2 (v.2.2.0)](https://daehwankimlab.github.io/hisat2/) (Kim et al. 2019) genome indexing will be applied. The reference genomes (.fasta) and the reference genome annotations (.gff) will be indexed with the hisat2-build command. The name of the genome indexes will be carefully chosen for easy downstream processes.

### Read mapping
[Hisat2 (v.2.2.0)](https://daehwankimlab.github.io/hisat2/) (Kim et al. 2019) conversion and sorting of alignment maps will be applied. Clean, well-trimmed fastq files will be mapped to their respective [reference genomes](https://www.ncbi.nlm.nih.gov/datasets/genome/). The input files will be the cleaned read files (paired .fastq) and the indexed genome. The input flags will be carefully chosen to produce the desired output. The output file expected from this command will be a sequence alignment map (SAM format). 
The SAM files will be converted to BAM files using the SAMtools view command. The BAM files will then be sorted using the [Samtools](https://www.htslib.org/doc/) sort command. SAMTools flag stat will be used to collect alignment statistics to ensure the quality of alignment before count matrices are generated. The input files for SAMtools (v.1.19) (Danecek et al. 2021) will be in SAM format while the output will be in sorted BAM format.

### Count Matrix
Count matrices will be generated using [StringTie (v.2.2.1)](https://ccb.jhu.edu/software/stringtie/) (Pertea et al. 2015) using the gene annotation files, reference genomes, and sorted BAM files for each sample. The output will be an assembly constructed from the RNAseq data compared to the reference genome and annotation files. Within StringTie, a prepD3.py3 python script (v.2.7.1) (Pertea et al. 2015) will be used to generate the count matrices from the gene transfer files for each sample. All the count matrices for all the samples will be joined together into one file.

### Differential Gene Expression (DGE)
Using the StringTie counts matrix, DESeq2 (v1.46.0.) (Love et al. 2014) will make an R object of the manually defined groups (treatments) and counts. Low-coverage genes will be filtered, and counts will be adjusted by normalization factors. An estimate of dispersion will be done. The DESeq2 output will be a new counts table (.csv), showing gene regulation

### Gene Enrichment and Molecular Network Analysis



