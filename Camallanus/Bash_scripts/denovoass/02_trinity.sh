#!/usr/bin/env bash

/usr/local/bin/Trinity \
  --seqType fq \
  --samples_file denovosamples.txt \
  --SS_lib_type FR \
  --max_memory 100G \
  --CPU 40 \
  --output ./trinity_out/
