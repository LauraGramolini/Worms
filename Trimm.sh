#!/usr/bin/env bash

###
#This is to use with parallel script, it works and produces 4 outputs (2 paired and 2 unpaired), and we should keep the .gz extension in the output
####


R1=$1
R2=$2
R3=$3
R4=$4

java -jar /home/laura/usr/share/java/trimmomatic-0.36.jar PE \
  -threads 10 \
  ./$R1 \
  ./$R2 \
  -baseout ./trim/$R3.fq.gz \
  ILLUMINACLIP:$TRIMMOMATIC_HOME/home/laura/usr/share/trimmomatic/TruSeq3-PE.fa:5:30:10 \
  SLIDINGWINDOW:5:5 \
  MINLEN:50 |&
  tee ./trim/$R3.log

######
#Then we use the trimmed for downstream
