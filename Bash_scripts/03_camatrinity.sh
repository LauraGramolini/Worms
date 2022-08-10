#!/usr/bin/env bash

###
#For denovo assembly I used Trinity, and 5 samples (1 for each of my condition) for a total of 41 Million reads
#I prepared a denovosamples.txt file that looks like this:
#Cond_A	Cond_A_rep1	LIB55_R1_1P.fastq.cor.fq.gz	LIB55_R1_2P.fastq.cor.fq.gz
#Cond_D	Cond_D_rep2	LIB41_R1_1P.fastq.cor.fq.gz	LIB41_R1_2P.fastq.cor.fq.gz
#Cond_E	Cond_E_rep2	LIB66_1P_7M.fastq.cor.fq.gz	LIB66_2P_7M.fastq.cor.fq.gz
#Cond_F	Cond_F_rep4	LIB51_R1_1P.fastq.cor.fq.gz	LIB51_R1_2P.fastq.cor.fq.gz
#Cond_G	Cond_G_rep1	LIB18_R1_1P.fastq.cor.fq.gz	LIB18_R1_2P.fastq.cor.fq.gz
#Using as inut the output from rcorrector.sh
#and then run Trinity like this:

/usr/local/bin/Trinity \
  --seqType fq \
  --samples_file denovosamples.txt \
  --SS_lib_type RF \
  --max_memory 10G \
  --CPU 10 \
  --output /SAN/Lauras_territory/Novaseq/cama/denovo/trinity_out/ #It creates the folder, just specify the name you want

##Then rename the headers to not have them mixed with other assemblies

sed 's/>TRINITY/>cama_TRINITY/g' Trinity.fasta > camadenovo.fasta
