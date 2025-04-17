# Toxin Identification Pipeline using Assembled and De Novo RNA Seq

[![]()]

[![]()]

This project focuses on the production of a pipeline to identify and (hopefully) dock proteins related to toxins by:

-
-

## Reference Based Pipeline:

### Required Packages 

|Package|Location|Use|
|  ------ | ------ | ------ |
|Conda|[Here](https://anaconda.org/anaconda/conda)| |
|BioConda|[Here]( )| |
|SRA-Tools|[Here]( )| |
|CutAdapt|[Here]( )| |
|Hisat2|[Here]( )| |
|Stringtie|[Here]( )| |
|Diamond|[Here]( )| |


```sh
conda create -n <env name> python=3.12
conda activate <env name>
conda deactivate
```

## De Novo Pipeline:

|Package|Location|Use|
|  ------ | ------ | ------ |
|Conda|[Here](https://anaconda.org/anaconda/conda)| |
|BioConda|[Here](https://bioconda.github.io/)| |
|SRA-Tools|[Here](https://github.com/ncbi/sra-tools)| |
|Trimmomatic|[Here](https://github.com/usadellab/Trimmomatic)| |
|Salmon|[Here](https://combine-lab.github.io/salmon/getting_started/)| |
|Trinity|[Here](https://combine-lab.github.io/salmon/getting_started/)| |

A majority of this project for De Novo was run on a server limited to 2-8 cores and 20-60GB DDR3 RAM 
This project can be run on minimum of these specs but was increased due to time constraints and limitations of those working on this pipeline.

