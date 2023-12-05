## First we find the ORF in our de novo assembly fasta file
TransDecoder.LongOrfs \
  -t ./SuperDuper.fasta &&

#######################The output is called longest_orfs.pep in the new directory and will be used for the hology finding##########################################


#with this now we match it for homology with Uniprot


blastp \
-query SuperDuper.fasta.transdecoder_dir/longest_orfs.pep  \
    -db /SAN/db/blastdb/uniprot/uniprot_sprot.fasta  -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 20 > ./annotation/blastp_trans.outfmt6 ; ###I'm not sure this will make them run in parallel to save time###

### Matches with Pfam

hmmscan \
  --cpu 20 \
  --domtblout ./annotation/pfam.domtblout \
  /SAN/db/interpro/interproscan-5.13-52.0/data/pfam/27.0/Pfam-A_newpress.hmm SuperDuper.fasta.transdecoder_dir/longest_orfs.pep  &&

###Then put all together for the prediction


/tools/TransDecoder-2.0.1/TransDecoder.Predict \
  -t SuperDuper.fasta \
  --retain_blastp_hits ./annotation/blastp_trans.outfmt6 \
  --retain_pfam_hits ./annotation/pfam.domtblout &&
##### From this you get 4 output and the .pep one is the important for the next steps (check the others, you can still do stuff with .bed)

###If it works, add this
#signalp -f short -n signalp.out camadenovo.fasta.transdecoder.pep
#And this
#tmhmm --short < camadenovo.fasta.transdecoder.pep > tmhmm.out


### Finally you can use eggnog to annotate
~/eggnog-mapper-2.1.7/emapper.py \
  --cpu 20 \
  -i SuperDuper.fasta.transdecoder.pep \
  -m diamond \
  -o ./eggnog \
  --excel
####It takes forever and the output is 3 files (I removed the log, I can put it back but it has to go at last!!!) : annotation, hits and see_orthologous that can be used in the next annotations. Change the base name of the output

