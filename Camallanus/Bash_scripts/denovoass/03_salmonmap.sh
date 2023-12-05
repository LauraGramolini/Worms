#!/usr/bin/env bash

for file in ./reads_xsalmon/*_mappedUnmapped.out.mate1.fq; do name=$(basename -a -s _mappedUnmapped.out.mate1.fq ${file});

salmon quant \
-i /SAN/Lauras_territory/cama_trans/salmonindexes \
--libType A \
--dumpEq \
--hardFilter \
--skipQuant \
-1 ./reads_xsalmon/${name}_mappedUnmapped.out.mate1.fq -2 ./reads_xsalmon/${name}_mappedUnmapped.out.mate2.fq \
-o salmon_out/${name}_out
done
