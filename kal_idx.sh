#!/bin/bash
../sw/kallisto_linux-v0.45.0/kallisto index -i Mus_musculus.GRCm38.cdna.all.fa.gz.idx Mus_musculus.GRCm38.cdna.all.fa.gz 

#zcat Mus_musculus.GRCm38.cdna.all.fa.gz | grep '>' | sed 's/gene_symbol:/\n/' | cut -d ' ' -f1 | paste - - | sed 's/>//' > genenames.tsv

zcat Mus_musculus.GRCm38.cdna.all.fa.gz | grep '>'  |   sed 's/gene:/\n/' | cut -d ' ' -f1 | paste - - | sed 's/>//' | cut -d '.' -f-2  > genenames.tsv 

