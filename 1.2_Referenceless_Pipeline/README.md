# De Novo Pipeline:
 
Although all of this poject was done in a Jupyter Server Environment (not using Conda), this process should work equally well within a standard conda environment. For this reason a YAML file was designed for those interested in running this through a simple conda environment and can be found [here](https://github.com/RIT-Environmental-Genomics/Toxin-Identification-Pipeline/blob/main/1.2_De_Novo_Pipeline/rnaseq_De_Novo.yml). Otherwise if you have previous experience with CONDA, be sure to download the packages from [here](https://github.com/RIT-Environmental-Genomics/Toxin-Identification-Pipeline/blob/main/README.md).

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


## 3. Trinity

Trinity was a backbone to forming this pipeline as a reference was needed to align 

## 4. Salmon
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

## 5: HISAT2

Top 50 TPM hits can be collected from the quant.sf output into a .txt format:

```sh
awk 'NR>1 {print $1, $5}' salmon_output/quant.sf | sort -k2,2nr | head -50 | cut -f1 -d' ' > top50_ids.txt
```

Then top 50 hits can be used with seqtk to subseq the origional Trinity.fasta to show only top hits:

```sh
seqtk subseq dir/Trinity.fasta dir/top50_ids.txt > Trinity50.fasta
```

HISAT2 can be used to aligntop hits back onto origional forward and reverse strands: 

Indexing:
```
hisat2-build top50_transcripts.fasta top50_index
```
Alignment
```sh
hisat2 -x top50_index \
  -1 forwardd_1.fastq -2 reverse_2.fastq \
  -S top50_aligned.sam \
  --threads 8
```

## 6: Convert and Filter

This step can now be skipped due to aligning using SEQTK and SALMON


## 7: Diamond [BlastX]

[Uniprot's](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz) Toxin database was used to identify genes related to toxins
(Due to the size of Uniprot's toxin database and the limitations of the hardware used up until this point [three laptops with less than 300 GB storage each], a server was used to contain and process the database into an unzipped fasta file)

``` sh
#unzip
gunzip uniprot_sprot.fasta.gz
```


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

Finished files were then compared with BLASTX from NCBI database to compare.






### Notes:

A majority of this project for De Novo was run on a server limited to 2-8 cores and 20-60GB DDR3 RAM 
This project can be run on minimum of these specs but was increased due to time constraints and limitations of those working on this pipeline.
