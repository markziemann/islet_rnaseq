#install.packages("devtools")
#library("devtools")
#devtools::install_github("markziemann/Mitch")
library("mitch")

# gene names
genenames<-read.table("../ref/mart_export_v96.txt",header=T,sep="\t")
gt<-unique(genenames[,c(1,3)])

# reactomes
genesets<-gmt_import("ReactomePathways_20180426.gmt")

dge1<-read.table("dge1_deseq.tsv")
dge2<-read.table("dge2_deseq.tsv")
dge3<-read.table("dge3_deseq.tsv")
dge4<-read.table("dge4_deseq.tsv")

w<-list("dge1"=dge1,"dge2"=dge2,"dge3"=dge3,"dge4"=dge4)

x<-mitch_import(w,DEtype="deseq2",geneTable=gt)

#run the analysis
res<-mitch_calc(x,genesets,resrows=250,bootstraps=1000,priority="effect")
mitch_plots(res,outfile="islet_reactome.pdf")
mitch_report(res,"islet_reactome.html")
