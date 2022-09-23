#!/usr/bin/env bash

###
#This comes after the Trimming, and it is also called by parallel.
#Here I correct the reads for illumina sequencing errors, it is the last step pf the general processing before assembling
##

R1=$1
R2=$2


perl ~/run_rcorrector.pl \
  -1 ./$R1 \
  -2 ./$R2 \
  -od ./corrected
