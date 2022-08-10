#!/usr/bin/env bash

## These are the annotation steps for our de novo assembly
#First we find the ORF in our de novo assembly fasta file


TransDecoder.LongOrfs \
  -t ./camadenovo.fasta &&

#######################The output is called longest_orfs.pep in the new directory and will be used for the holmogy finding##########################################


#with this now we match it for homology with Uniprot ...

blastp \
-query camadenovo.fasta.transdecoder_dir/longest_orfs.pep  \
    -db /SAN/db/blastdb/uniprot/uniprot_sprot.fasta  -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 20 > ./blastp_trans.outfmt6  ;

### ... with Pfam ...

hmmscan \
  --cpu 8 \
  --domtblout pfam.domtblout \
  /SAN/db/interpro/interproscan-5.13-52.0/data/pfam/27.0/Pfam-A_newpress.hmm camadenovo.fasta.transdecoder_dir/longest_orfs.pep ;


###Then put all together for the prediction

/tools/TransDecoder-2.0.1/TransDecoder.Predict \
  -t camadenovo.fasta \
  --retain_blastp_hits blastp_trans.outfmt6 \
  --retain_pfam_hits pfam.domtblout ;

##### From this you get 4 output and the camadenovo.fasta.transdecoder.pep is the important for the next step (check the others, you can still do stuff with .bed), the gff3 can also
##be included in the eggnog step

##################################################As additional prediction we could use signalp and tmhmm ##############################################
### ...with signal peptide prediction ...
signalp -f short -n signalp.out camadenovo.fasta.transdecoder.pep ;

### ...with transmembrane prediction ...
tmhmm --short < camadenovo.fasta.transdecoder.pep > tmhmm.out ;
#########################################################################################################################################################

### Finally you can use eggnog to annotate, starting from the Transdecoder predict output
~/eggnog-mapper-2.1.7/emapper.py \
  --cpu 4 \
  -i camadenovo.fasta.transdecoder.pep \
  -m diamond \
  -o ./eggnog \
  --decorate_gff camadenovo.fasta.transdecoder.gff3 \
  --excel

####The last part is rerunning to be sure,
#Outuputs files: _ annotation, _hits and _see_orthologous that can be used in the next annotations

