## De Novo Pipeline:

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


Although all of this poject was done in a Jupyter Server Environment (not using Conda), this process should work equally well within a standard conda environment. For this reason a YAML file was designed for those interested in running this through a simple conda environment and can be found [here](https://github.com/RIT-Environmental-Genomics/Toxin-Identification-Pipeline/blob/main/1.2_De_Novo_Pipeline/rnaseq_De_Novo.yml).

|Package|Use|
| ------ | ------ |
|[SRA-Tools](https://github.com/ncbi/sra-tools)| |
|[Trimmomatic](https://github.com/usadellab/Trimmomatic)| |
|[Salmon](https://combine-lab.github.io/salmon/getting_started/)| |
|[Trinity](https://combine-lab.github.io/salmon/getting_started/)| |




A majority of this project for De Novo was run on a server limited to 2-8 cores and 20-60GB DDR3 RAM 
This project can be run on minimum of these specs but was increased due to time constraints and limitations of those working on this pipeline.
