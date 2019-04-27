#!/bin/bash

#for FQZ in *fastq.gz ; do
# ../sw/skewer -q 10 $FQZ -t 16
#done

for FQ in *fastq ; do
  ../sw/kallisto_linux-v0.45.0/kallisto quant \
  -i ../ref/Mus_musculus.GRCm38.cdna.all.fa.gz.idx \
  -o ${FQ}_kal -t 16 --single -l 150 -s 10 $FQ
done

for TSV in */*abundance.tsv ; do
  NAME=$(echo $TSV | cut -d '_' -f1) ; cut -f1,4 $TSV | sed 1d | sed "s/^/${NAME}\t/"
done > 3col.tsv
