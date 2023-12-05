#!/usr/bin/env bash


for file in /SAN/Lauras_territory/Novaseq/cama/raw_reads/*_R1.fastq.gz; do name=$(basename -a -s _R1.fastq.gz ${file});

STAR \
--runThreadN 40 \
--genomeDir ./starindex_super \
--readFilesIn /SAN/Lauras_territory/Novaseq/cama/raw_reads/${name}_R1.fastq.gz /SAN/Lauras_territory/Novaseq/cama/raw_reads/${name}_R2.fastq.gz \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix ./star_out_superloop/${name}_mapped
done
