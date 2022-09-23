#!/usr/bin/env bash

###After rcorrector I used 1 sample only (each copepod sample is made of 8-10 copepods) for 400 M reads

/usr/local/bin/Trinity \
  --seqType fq \
  --left LIB27_R1_1P.fastq.cor.fq.gz \
  --right LIB27_R1_2P.fastq.cor.fq.gz \
  --SS_lib_type RF \
  --max_memory 200G \
  --CPU 40 \
  --output ./trinity_outLIB27

##Then rename the headers to not have them mixed with other assemblies

sed 's/>TRINITY/>cop_TRINITY/g' Trinity.fasta > copepoddenovo.fasta
