#!/usr/bin/env bash

###
#I've used trimmomatic to trim my raw reads as preprocessing before assemblying the de novo transcriptomes
#I've used this script for Cmallanus and the copepod raw reads
#This script is called by parallel so my samples are processed all together. It produces 4 outputs (2 paired and 2 unpaired), and we need the Paired pair for downstream
####


R1=$1
R2=$2
R3=$3

java -jar /home/laura/usr/share/java/trimmomatic-0.36.jar PE \
  -threads 10 \
  ./$R1 \
  ./$R2 \
  -baseout ./trim/$R3.fq.gz \
  ILLUMINACLIP:$TRIMMOMATIC_HOME/home/laura/usr/share/trimmomatic/TruSeq3-PE.fa:5:30:10 \
  SLIDINGWINDOW:5:5 \
  MINLEN:50 |&
  tee ./trim/$R3.log
