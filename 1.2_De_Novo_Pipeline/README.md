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
salmon index -t trinity_out_dir/Trinity.fasta \
 -i salmon/trinity_index
```

```sh
salmon quant -i trinity_index \
  -l A \
  -1  dir/forward_1.fastq \
  -2  dir/reverse_2.fastq \
  -p 8 \
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


## 4: Convert and Filter

Sam files were then converted to BAM files using SAMTools:

```sh
samtools view -S -b input.sam > output.bam
```

or
```
samtools sort -o sorted_output.bam output.bam
```

Stringtie was then used to convert BAM files into GTF files:

```sh
stringtie input.bam -o output.gtf
```
Finally GTF files were converted into GFF3 using GFFREAD:
```sh
gffread input.gtf -T -o- | gffread - -E -o output.gff3
```

GFF3 files were then filtered based on Transcripts Per Million ```TPM``` with only top 25 percentile of most abundant being left after filtering was completed

After filtering, [Uniprot's](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz) Toxin database was used to identify genes related to toxins
(Due to the size of Uniprot's toxin database and the limitations of the hardware used up until this point [three laptops with less than 300 GB storage each], a server was used to contain and process the database into an unzipped fasta file)

``` sh
#unzip
gunzip uniprot_sprot.fasta.gz
```

## 5: Diamond [BlastX]
create diamond database from uniprot fasta from only proteins from the species 

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
	-s title
```

Finished files were then compared with NCBI database to confirm accuracy of identified genes.






### Notes:

A majority of this project for De Novo was run on a server limited to 2-8 cores and 20-60GB DDR3 RAM 
This project can be run on minimum of these specs but was increased due to time constraints and limitations of those working on this pipeline.
