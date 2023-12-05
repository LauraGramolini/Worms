#!/usr/bin/env bash


for file in ./*R1.fastq.gz; do name=$(basename -a -s _R1.fastq.gz ${file});

STAR \
--runThreadN 40 \
--genomeDir /SAN/Lauras_territory/dual_ref/starindex_superdual/ \
--readFilesIn ./${name}_R1.fastq.gz ./${name}_R2.fastq.gz \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix ./star_superdual/${name}_mapped
done
