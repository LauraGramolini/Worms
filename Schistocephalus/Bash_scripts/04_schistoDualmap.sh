#!/usr/bin/env bash

##For schistocephalus-copepod dual-seq samples, I had the assembly (03_coptrinity.sh) and I have to concatenate it with the schistocephalus genome

cat /SAN/Lauras_territory/Novaseq/copepod/assembly/LIB27/trinity_outLIB27.Trinity.fasta /SAN/Lauras_territory/schisto_genome/schistocephalus_solidus.PRJEB527.WBPS16.genomic.fa >> LIB27schistocop_ref.fasta

##This is to create the indexes for mapping for this concatenated ref (copepod transcriptme + Schistocephalus genome)

STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir ./LIB27star_index \
--genomeFastaFiles /SAN/Lauras_territory/dual_ref/LIB27schistocop_ref.fasta \
--limitGenomeGenerateRAM 310000000000


# 195GB=210000000000bytes

###Then I use the indexes to map the dual-seq samples

for file in PATH/*R1.fastq.gz; do name=$(basename -a -s _R1.fastq.gz ${file});

STAR \
--runThreadN 10 \
--genomeDir PATH/dual_ref/LIB27star_index \
--readFilesIn ./${name}_R1.fastq.gz ./${name}_R2.fastq.gz \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--sjdbGTFfile /SAN/Lauras_territory/schisto_genome/schisto_annotation.gtf \
--outFileNamePrefix ./star_out/${name}_mapped
done

##Then I move the .bam files where the other .bam from single sequencing schistocephalus data are

mv *.bam PATH/AllschistoBam

##Then I run FeatureCounts to count the reads and get the count matrix
#Scripts in the R folder
