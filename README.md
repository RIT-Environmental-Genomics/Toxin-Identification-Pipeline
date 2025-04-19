# Toxin Identification Pipeline using Assembled and De Novo RNA Seq

[![]()]

[![]()]

This project focuses on the production of a pipeline to identify and (hopefully) dock proteins related to toxins by:

-
-

## Reference Based Pipeline: 

### Required Packages 

|Repository| 
|  ------ | 
|[Conda](https://anaconda.org/anaconda/conda)| 
|[Conda-Forge](https://conda-forge.org/)| 
|[BioConda](https://bioconda.github.io/)| 



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

### Step 1: Fastqc

SRA-Tools prefetch command was used to download RNA seq in SRA format
```sh
prefetch <Assencion ID>
```

```sh
Fast-dump --split-files <Assencion ID>
```





## De Novo Pipeline:

|  ------ |  
|[Conda](https://anaconda.org/anaconda/conda)| 
|[Conda-Forge](https://conda-forge.org/)| 
|[BioConda](https://bioconda.github.io/)| 

|Package|Use|
|[SRA-Tools](https://github.com/ncbi/sra-tools)| |
|[Trimmomatic](https://github.com/usadellab/Trimmomatic)| |
|[Salmon](https://combine-lab.github.io/salmon/getting_started/)| |
|[Trinity](https://combine-lab.github.io/salmon/getting_started/)| |

A majority of this project for De Novo was run on a server limited to 2-8 cores and 20-60GB DDR3 RAM 
This project can be run on minimum of these specs but was increased due to time constraints and limitations of those working on this pipeline.

