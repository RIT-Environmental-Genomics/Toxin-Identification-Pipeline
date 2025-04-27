<p>
  
<h1>Toxin Identification Pipeline of RNA Seq data using Reference and <br >
De Novo genome sequences </h1>



![Species](other/pngs/Species.png)

![Pipeline](https://github.com/RIT-Environmental-Genomics/Toxin-Identification-Pipeline/blob/main/other/pngs/Pipeline.png)

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

<br >

### Pipelines:
|Pipeline|
| ------ |
|[Reference](https://github.com/RIT-Environmental-Genomics/Toxicology/blob/main/1.1_Reference_Pipeline/)|
|[De Novo](https://github.com/RIT-Environmental-Genomics/Toxicology/tree/main/1.2_De_Novo_Pipeline)|

<br >

### Conda Environment Repositories: 
  
|Repository| 
|  ------ | 
|[Conda](https://anaconda.org/anaconda/conda)| 
|[Conda-Forge](https://conda-forge.org/)| 
|[BioConda](https://bioconda.github.io/)| 

<br >

## Package

### All Pipelines
|Package|Use|
|  ------ | ------ |
|[SRA-Tools](https://github.com/ncbi/sra-tools)| |
|[Fastqc](https://bioconda.github.io/recipes/fastqc/README.html)||
|[CutAdapt](https://anaconda.org/bioconda/cutadapt)| |
|[Hisat2](https://anaconda.org/bioconda/hisat2)| |
|[Stringtie](https://anaconda.org/bioconda/stringtie)| |
|[Samtools](https://anaconda.org/bioconda/samtools)| |
|[Gffread](https://anaconda.org/bioconda/gffread)| |
|[Seqtk](https://anaconda.org/bioconda/seqtk)| |
|[Diamond](https://anaconda.org/bioconda/diamond)| |  

### De Novo Pipeline Only
|Package|Use|
|  ------ | ------ |
|[Trimmomatic](https://github.com/usadellab/Trimmomatic)| |
|[Salmon](https://combine-lab.github.io/salmon/getting_started/)| |
|[Trinity](https://combine-lab.github.io/salmon/getting_started/)| |


|Results|
|:-:|
|[FASTQC](https://rit-environmental-genomics.github.io/Toxin-Identification-Pipeline/Results/FASTQC/index.html)|


</p>
