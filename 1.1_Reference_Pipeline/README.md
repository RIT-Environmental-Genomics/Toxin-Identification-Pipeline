
### Required Packages 

|Repository| 
|  ------ | 
|[Conda](https://anaconda.org/anaconda/conda)| 
|[Conda-Forge](https://conda-forge.org/)| 
|[BioConda](https://bioconda.github.io/)| 

## 1: Create new Conda Environment

Through Conda environment, bioconda and conda-forge repositories were added:

```sh
add bioconda and conda-forge repositories
conda config –add channels bioconda
conda config –add channels conda-forge
```

A new Conda environment can then be used in any of the three following ways

```sh
conda create -n <env name> python sra-tools cutadapt hisat2 stringtie samtools bedtools seqtk diamond
conda activate <env name>
conda deactivate
```
or 
```sh
conda create -f test.yml
```
or
```sh
conda create -n <env name>
conda activate <env name>
conda install <package>
conda deactivate
```



|Package|Use|
|  ------ | ------ |
|[SRA-Tools](https://github.com/ncbi/sra-tools)| |
|[Fastqc]()||
|[CutAdapt]()| |
|[Hisat2]()| |
|[Stringtie]()| |
|[Samtools]()| |
|[Bedtools]()| |
|[Seqtk]()| |
|[Diamond]()| |

## 2: Prefetch and Fastqc

SRA-Tools prefetch command was used to download RNA seq in SRA format
```sh
prefetch <Assencion ID>
```
Fast-dump was then used to split the file into forward and reverse strands and FASTQ program's fasatqc command was used to analyze quality of FASTA files
```sh
Fast-dump --split-files <SRR ID>

fastqc <SRRID_1.fasta>
fastqc <SRRID_2.fasta>
```
Forward and reverse strands are identified using _1 (Forward) and _2 (Reverse). 

## 3: Trimmomatic or CutAdapt
CutAdapt was used to cut any remaining adapters or unreliable sequences found in FASTQC, however Trimmomatic was used in place of CutAdapt for use in Trinity

```sh
cutadapt -a <Sequence in ATCG format> -o output.fasta input.fasta
```

## 4: FASTA Alignment
Hisat2 was then used to build a genome-index based on a reference genome and aligned to forward and reverse reads to provide a " schaffold " of the fully alligned strands in SAM format
```sh
hisat2-build -p <threads> <Reference_Genome>_genomic.fna genome-index

hisat2 -p <threads> -x genome-index --dta -1 <Forward Strand>.fastq -2 <Reverse_Strand>.fastq -S <aligned_strands>.sam
```

Sam files were then converted to BAM files

Bam files converted to GTF files

GTF files converted to GFF3 files

GFF3 files were then filtered based on Transcripts Per Million ```TPM``` with only top 25 percentile of most abundant being left after filtering was completed

After filtering, [Uniprot's](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz) Toxin database was used to 
``` sh
#unzip
gunzip uniprot_sprot.fasta.gz

#


```
## Diamond [BlastX]
create diamond database from uniport fasta from only proteins from the species 
```sh
grep -A 1 ">.*<species taxa>" uniprot_sprot.fasta > <species>.fasta 
diamond makedb --in <species>.fasta -d <species>_db
```

Example:
```sh
grep -A 1 ">.*Pseudonaja textilis" uniprot_sprot.fasta > snake_only_sprot.fasta 
diamond makedb --in snake_only_sprot.fasta -d snake_sprot_db
```
Diamond BlastX was then done using merged FASTA files and _db.dmnd file
```sh
diamond blastx \
	-d <database>.dmnd \
	-q <filtered_merged_sequences>.fasta \
	-o <Species_Results>.m8 \
	-f 6 qseqid sseqid pident evalue bitscore stitle\
	-p <processors/threads> \
#if you need Sensitive:
	--sensitive \
	-k 1 <top searches> \
#if you need to write to a temporary folder
	-t /tmp/diamond_temp 
```

