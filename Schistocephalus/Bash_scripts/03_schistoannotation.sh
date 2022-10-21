###We extract the transcripts from the ref genome of Schisto:

gffread -F -w schisto_transcripts.fasta -g schisto_genome.fasta schistocephalus_solidus.PRJEB527.WBPS16.annotations.gff3

###This produced 2 files: schisto_genome.fasta.fai (with indexes) and schisto_transcripts.fasta
#Now we can use the second with TransDecoder

TransDecoder.LongOrfs -t schisto_transcripts.fasta

###It produces the folder "schisto_transcripts.fasta.transdecoder_dir/" with the l-orf file
#Searching for omologies with blast


blastp \
-query schisto_transcripts.fasta.transdecoder_dir/longest_orfs.pep  \
    -db /SAN/db/blastdb/uniprot/uniprot_sprot.fasta  -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 10 > ./schisto_blastp.outfmt6  \

### Matches with Pfam

hmmscan \
  --cpu 8 \
  --domtblout schisto_pfam.domtblout \
  /SAN/db/interpro/interproscan-5.13-52.0/data/pfam/27.0/Pfam-A_newpress.hmm schisto_transcripts.fasta.transdecoder_dir/longest_orfs.pep \

###Then put all together for the prediction

/tools/TransDecoder-2.0.1/TransDecoder.Predict \
  -t schisto_transcripts.fasta \
  --retain_blastp_hits schisto_blastp.outfmt6 \
  --retain_pfam_hits schisto_pfam.domtblout \

### From this you get 4 output and the .pep one is the important for the next steps (check the others, you can still do stuff with .bed)


### Finally you can use eggnog to annotate
~/eggnog-mapper-2.1.7/emapper.py \
  --cpu 4 \
  -i schisto_genome.fasta.transdecoder.pep \
  -m diamond \
  -o ./eggnog |& tee ./eggnog.log

##We obtain our annotation file: "eggnong.emapper.annotations" with the GO terms
