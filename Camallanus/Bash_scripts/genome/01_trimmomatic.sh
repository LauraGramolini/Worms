#!/usr/bin/env bash

java -jar /home/laura/usr/share/java/trimmomatic-0.36.jar PE \
  -threads 10 \
  /SAN/Lauras_territory/cama_genome/Miseq_2023/Run1/S0001-L3_S1_L001_R1_001.fastq.gz \
  /SAN/Lauras_territory/cama_genome/Miseq_2023/Run1/S0001-L3_S1_L001_R2_001.fastq.gz \
  -baseout /SAN/Lauras_territory/cama_genome/Miseq_2023/trimrun1/camagentrim.fq.gz \
  ILLUMINACLIP:$TRIMMOMATIC_HOME/home/laura/usr/share/trimmomatic/TruSeq3-PE.fa:5:30:10 \
  SLIDINGWINDOW:4:15 \
  MINLEN:50 \
  -phred33 |&
  tee /SAN/Lauras_territory/cama_genome/Miseq_2023/trimrun1/camagentrim.log
