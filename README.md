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
|[BioConda](https://bioconda.github.io/)| 




```sh
add bioconda and conda-forge repositories
conda config –add channels bioconda
conda config –add channels conda-forge
```


```sh
conda create -n <env name> python=3.12
conda activate <env name>
conda deactivate
```

|Package|Use|
|  ------ | ------ |
|[SRA-Tools](https://github.com/ncbi/sra-tools)| |
|[CutAdapt]()| |
|[Hisat2]()| |
|[Stringtie]()| |
|[Samtools]()| |
|[Bedtools]()| |
|[Seqtk]()| |
|[Diamond]()| |

## De Novo Pipeline:

|Package|Use|
|  ------ | ------ |
|[Conda](https://anaconda.org/anaconda/conda)| |
|[BioConda](https://bioconda.github.io/)| |
|[SRA-Tools](https://github.com/ncbi/sra-tools)| |
|[Trimmomatic](https://github.com/usadellab/Trimmomatic)| |
|[Salmon](https://combine-lab.github.io/salmon/getting_started/)| |
|[Trinity](https://combine-lab.github.io/salmon/getting_started/)| |

A majority of this project for De Novo was run on a server limited to 2-8 cores and 20-60GB DDR3 RAM 
This project can be run on minimum of these specs but was increased due to time constraints and limitations of those working on this pipeline.

