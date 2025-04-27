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


## 2. Trinity

## 3. Salmon
Instead of using Trinity's built in Salmon quantification, Salmon was used at the final steps manually to quantify data as follows:

Trinity.fasta was indexed (just as in HISAT2) 
```
salmon index -t Trinity.fasta -i trinity_index
```

Forward and Reverse strands from SRA file were then used to quantify data:

```sh
salmon quant -i trinity_index -l A \
             -1 forward_1.fastq -2 reverse_2.fastq \
             -p 8\
             -o salmon_output
```

Finally, TPM was extracted to help filter for remaining steps

```sh
abundance_estimates_to_matrix.pl \
    --est_method salmon \
    --gene_trans_map Trinity.fasta.gene_trans_map \
    --name_sample_by_basedir \
    salmon_output/quant.sf
```

### 

A majority of this project for De Novo was run on a server limited to 2-8 cores and 20-60GB DDR3 RAM 
This project can be run on minimum of these specs but was increased due to time constraints and limitations of those working on this pipeline.
