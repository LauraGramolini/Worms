#!/bin/bash

##For Schistocephalus pure samples I directly mapped against the ref genome from Wormbase parasite
#This is to create the indexes previous mapping

STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir ./indexes \
--genomeFastaFiles /SAN/Lauras_territory/schisto_genome/schistocephalus_solidus.PRJEB527.WBPS16.genomic.fa \
--sjdbGTFfile /SAN/Lauras_territory/schisto_genome/schistocephalus_solidus.PRJEB527.WBPS16.annotations.gff3 \
--sjdbGTFtagExonParentTranscript Parent \
--limitGenomeGenerateRAM 210000000000


# 195GB=210000000000bytes

###Then the indexes were used to map the raw reads from schistocephalus pure samples

for file in PATH/*R1.fastq.gz; do name=$(basename -a -s _R1.fastq.gz ${file});

STAR \
--runThreadN 10 \
--genomeDir /SAN/Lauras_territory/schisto_genome/star_indexes/ \
--readFilesIn ./${name}_R1.fastq.gz ./${name}_R2.fastq.gz \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--sjdbGTFfile /SAN/Lauras_territory/schisto_genome/schisto_annotation.gtf \
--quantMode GeneCounts \
--outFileNamePrefix ./star_out/${name}_mapped
done


###Then all the .bam files were moved to a common folder together with the .bam from dual-seq

mv PATH/*.bam PATH/AllschistoBams/
