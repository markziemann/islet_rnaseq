# islet_rnaseq

In this project, the scripts were run in the following order:

1. kal_idx.sh to index the mouse transcriptome and output a table of gene accessions/symbols

2. run_kal.sh to perform quality trimming and map the reads to the mouse transcriptome

3. deseq2.R to import expression counts into R and perform differential analysis.

4. run_mitch.R to perform enrichment analysis
