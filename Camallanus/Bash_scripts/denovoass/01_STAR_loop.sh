#!/usr/bin/env bash


for file in ../reads/corrected/*_R1_1P.fastq.cor.fq.gz; do name=$(basename -a -s _R1_1P.fastq.cor.fq.gz ${file});

STAR \
--runThreadN 10 \
--genomeDir ./indexes_fish \
--readFilesIn ../reads/corrected/${name}_R1_1P.fastq.cor.fq.gz ../reads/corrected/${name}_R1_2P.fastq.cor.fq.gz \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--sjdbGTFfile ./fish_annot.gff \
--limitSjdbInsertNsj 2428996 \
--outFileNamePrefix ./star_out_loop/${name}_mapped
done
