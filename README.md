# Toxin Identification Pipeline of RNA Seq data using Reference and De Novo genome sequences

[![]()]

[![]()]

This project focuses on the production of a pipeline to identify and (hopefully) dock proteins related to toxins by:

- Uploading RNASeq Data from NCBI
- Analyzing sequence quality and trimming both adapters and unreliable end sequences
- Convert sequences into SAM using:
  - Scaffolding them to a reference genome
  - Assembled as a De Novo sequence using Trinity Assembler
-  Convert SAM into GFF3 and filter out all but top 25% sequence abundance
-  Filter forward and reverse Sequences using GFF3
-  BLASTX sequence to Toxin Database [Uniprot](https://www.uniprot.org/) based on species 
-  Compare BLASTX results to 

|Repository| 
|  ------ | 
|[Conda](https://anaconda.org/anaconda/conda)| 
|[Conda-Forge](https://conda-forge.org/)| 
|[BioConda](https://bioconda.github.io/)| 


### All Pipelines
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

### De Novo Pipeline Only
|Package|Use|
|  ------ | ------ |
|[Trimmomatic](https://github.com/usadellab/Trimmomatic)| |
|[Salmon](https://combine-lab.github.io/salmon/getting_started/)| |
|[Trinity](https://combine-lab.github.io/salmon/getting_started/)| |
