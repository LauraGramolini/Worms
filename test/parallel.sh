#!/usr/bin/env bash

###
#This works. It calls the sortmenra script for all the samples in parallel. But the output is only one sample. I guess the problem is the sortmerna.sh
#For now is in pause
###
parallel  --xapply ./Trimm.sh {1} {2} {1.} {2.} ::: *R1_001.fastq ::: *R2_001.fastq
