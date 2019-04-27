library("tidyverse")
library("reshape2")
library("DESeq2")
library("gplots")

# read counts
tmp<-read.table("3col.tsv",header=F)
x<-as.matrix(acast(tmp, V2~V1, value.var="V3", fun.aggregate = sum))

# get gene names
gt<-read.table("../ref/genenames.tsv")
xx<-merge(x,gt,by.x=0,by.y="V1")
xx$Row.names=NULL
xxx<-aggregate(. ~ V2,xx,sum)
rownames(xxx)<-xxx$V2
xxx$V2=NULL
x<-round(xxx)

# ortholog file for mapping reactomes
orth<-read.table("../ref/mart_export_v96.txt",sep="\t",header=T)
orth$Transcript.stable.ID=NULL

# samplesheet curation
samples<-as.data.frame(colnames(x))
colnames(samples)<-"sample"
samples$hgww<-as.numeric(grepl("1$",samples$sample))
samples$ngww<-as.numeric(grepl("2$",samples$sample))
samples$hgmm<-as.numeric(grepl("3$",samples$sample))
samples$ngmm<-as.numeric(grepl("4$",samples$sample))
rownames(samples)<-samples$sample
samples$sample=NULL
samples$batch<-as.factor(sapply(strsplit(rownames(samples),"-"),"[[",2))

####################################################
# contrast ngww vs hgww
####################################################
des<-subset( samples, ngww == 1 | hgww ==1 )
y<-x[,which(colnames(x) %in% rownames(des))]
y<-y[which(rowSums(y)/ncol(y)>=(10)),]
dds <- DESeqDataSetFromMatrix(countData = y , colData = des, design = ~ batch + hgww )
res <- DESeq(dds)
z<- results(res)
vsd <- vst(dds, blind=FALSE)
zz<-cbind(as.data.frame(z),assay(vsd))
dge<-as.data.frame(zz[order(zz$pvalue),])
write.table(dge,file="dge1_deseq.tsv",quote=F,sep="\t")
rnk<-as.data.frame( sign(dge$log2FoldChange) * (-log(dge$pvalue + 1E-307) ))
rownames(rnk)<-row.names(dge)
colnames(rnk)="Score"
rnk<-unique(merge(rnk,orth,by.x=0,by.y="Gene.stable.ID"))
rnk$Row.names=NULL
rnk<-aggregate(. ~ Human.gene.name,rnk,sum)
rnk<-rnk[-which(rnk$Human.gene.name==""),]
rownames(rnk)<-rnk$Human.gene.name
rnk$Human.gene.name=NULL
write.table(rnk,file="dge1.rnk",sep='\t',quote=F)
rnk1<-rnk 
dge1<-dge

#some plots
pdf("dge1_plots.pdf")
sig<-subset(dge,padj<0.05)
SIG=nrow(sig)
DN=nrow(subset(sig,log2FoldChange<0))
UP=nrow(subset(sig,log2FoldChange>0))
HEADER=paste("NGWW vs HGWW:", SIG , "DGEs,", UP ,"upregulated,", DN, "downregulated")
plot(log2(dge$baseMean),dge$log2FoldChange,cex=0.6, xlab="log2 base mean", ylab="log2 fold change" ,pch=19,col="#838383")
points(log2(sig$baseMean),sig$log2FoldChange,cex=0.6,pch=19,col="red")
mtext(HEADER)
top<-head(sig,20)
#text(log2(top$baseMean)+1, top$log2FoldChange, labels = rownames(top),cex=0.7)
#volcano plot
plot(dge$log2FoldChange, -log2(dge$pvalue) ,cex=0.6, xlim=c(-4,6),xlab="log2 fold change", ylab="-log2 p-value" ,pch=19,col="#838383")
points(sig$log2FoldChange, -log2(sig$pvalue),cex=0.6,pch=19,col="red")
#text(top$log2FoldChange+0.5, -log2(top$pvalue), labels = rownames(top),cex=0.7)
mtext(HEADER)
# top N gene heatmap
colfunc <- colorRampPalette(c("blue", "white", "red"))
heatmap.2(  as.matrix(dge[1:100,c(7:ncol(dge))]), col=colfunc(25),scale="row",
 trace="none",margins = c(6,6), cexRow=.4, main="Top 100 genes")
dev.off()

####################################################
# contrast ngmm vs hgmm
####################################################
des<-subset( samples, ngmm == 1 | hgmm ==1 )
y<-x[,which(colnames(x) %in% rownames(des))]
y<-y[which(rowSums(y)/ncol(y)>=(10)),]
dds <- DESeqDataSetFromMatrix(countData = y , colData = des, design = ~ batch + hgmm )
res <- DESeq(dds)
z<- results(res)
vsd <- vst(dds, blind=FALSE)
zz<-cbind(as.data.frame(z),assay(vsd))
dge<-as.data.frame(zz[order(zz$pvalue),])
write.table(dge,file="dge2_deseq.tsv",quote=F,sep="\t")
rnk<-as.data.frame( sign(dge$log2FoldChange) * (-log(dge$pvalue + 1E-307) ))
rownames(rnk)<-row.names(dge)
colnames(rnk)="Score"
rnk<-unique(merge(rnk,orth,by.x=0,by.y="Gene.stable.ID"))
rnk$Row.names=NULL
rnk<-aggregate(. ~ Human.gene.name,rnk,sum)
rnk<-rnk[-which(rnk$Human.gene.name==""),]
rownames(rnk)<-rnk$Human.gene.name
rnk$Human.gene.name=NULL
write.table(rnk,file="dge2.rnk",sep='\t',quote=F)
rnk2<-rnk
dge2<-dge

#some plots
pdf("dge2_plots.pdf")
sig<-subset(dge,padj<0.05)
SIG=nrow(sig)
DN=nrow(subset(sig,log2FoldChange<0))
UP=nrow(subset(sig,log2FoldChange>0))
HEADER=paste("NGMM vs HGMM:", SIG , "DGEs,", UP ,"upregulated,", DN, "downregulated")
plot(log2(dge$baseMean),dge$log2FoldChange,cex=0.6, xlab="log2 base mean", ylab="log2 fold change" ,pch=19,col="#838383")
points(log2(sig$baseMean),sig$log2FoldChange,cex=0.6,pch=19,col="red")
mtext(HEADER)
top<-head(sig,20)
#text(log2(top$baseMean)+1, top$log2FoldChange, labels = rownames(top),cex=0.7)
#volcano plot
plot(dge$log2FoldChange, -log2(dge$pvalue) ,cex=0.6, xlim=c(-4,6),xlab="log2 fold change", ylab="-log2 p-value" ,pch=19,col="#838383")
points(sig$log2FoldChange, -log2(sig$pvalue),cex=0.6,pch=19,col="red")
#text(top$log2FoldChange+0.5, -log2(top$pvalue), labels = rownames(top),cex=0.7)
mtext(HEADER)
# top N gene heatmap
colfunc <- colorRampPalette(c("blue", "white", "red"))
heatmap.2(  as.matrix(dge[1:100,c(7:ncol(dge))]), col=colfunc(25),scale="row",
 trace="none",margins = c(6,6), cexRow=.4, main="Top 100 genes")
dev.off()

####################################################
# contrast ngmm vs ngww
####################################################
des<-subset( samples, ngmm == 1 | ngww ==1 )
y<-x[,which(colnames(x) %in% rownames(des))]
y<-y[which(rowSums(y)/ncol(y)>=(10)),]
dds <- DESeqDataSetFromMatrix(countData = y , colData = des, design = ~ batch + ngww )
res <- DESeq(dds)
z<- results(res)
vsd <- vst(dds, blind=FALSE)
zz<-cbind(as.data.frame(z),assay(vsd))
dge<-as.data.frame(zz[order(zz$pvalue),])
write.table(dge,file="dge3_deseq.tsv",quote=F,sep="\t")
rnk<-as.data.frame( sign(dge$log2FoldChange) * (-log(dge$pvalue + 1E-307) ))
rownames(rnk)<-row.names(dge)
colnames(rnk)="Score"
rnk<-unique(merge(rnk,orth,by.x=0,by.y="Gene.stable.ID"))
rnk$Row.names=NULL
rnk<-aggregate(. ~ Human.gene.name,rnk,sum)
rnk<-rnk[-which(rnk$Human.gene.name==""),]
rownames(rnk)<-rnk$Human.gene.name
rnk$Human.gene.name=NULL
write.table(rnk,file="dge3.rnk",sep='\t',quote=F)
rnk3<-rnk
dge3<-dge

#some plots
pdf("dge3_plots.pdf")
sig<-subset(dge,padj<0.05)
SIG=nrow(sig)
DN=nrow(subset(sig,log2FoldChange<0))
UP=nrow(subset(sig,log2FoldChange>0))
HEADER=paste("NGMM vs NGWW:", SIG , "DGEs,", UP ,"upregulated,", DN, "downregulated")
plot(log2(dge$baseMean),dge$log2FoldChange,cex=0.6, xlab="log2 base mean", ylab="log2 fold change" ,pch=19,col="#838383")
points(log2(sig$baseMean),sig$log2FoldChange,cex=0.6,pch=19,col="red")
mtext(HEADER)
top<-head(sig,20)
#text(log2(top$baseMean)+1, top$log2FoldChange, labels = rownames(top),cex=0.7)
#volcano plot
plot(dge$log2FoldChange, -log2(dge$pvalue) ,cex=0.6, xlim=c(-4,6),xlab="log2 fold change", ylab="-log2 p-value" ,pch=19,col="#838383")
points(sig$log2FoldChange, -log2(sig$pvalue),cex=0.6,pch=19,col="red")
#text(top$log2FoldChange+0.5, -log2(top$pvalue), labels = rownames(top),cex=0.7)
mtext(HEADER)
# top N gene heatmap
colfunc <- colorRampPalette(c("blue", "white", "red"))
heatmap.2(  as.matrix(dge[1:100,c(7:ncol(dge))]), col=colfunc(25),scale="row",
 trace="none",margins = c(6,6), cexRow=.4, main="Top 100 genes")
dev.off()

####################################################
# contrast hgmm vs hgww
####################################################
des<-subset( samples, hgmm == 1 | hgww ==1 )
y<-x[,which(colnames(x) %in% rownames(des))]
y<-y[which(rowSums(y)/ncol(y)>=(10)),]
dds <- DESeqDataSetFromMatrix(countData = y , colData = des, design = ~ batch + hgww )
res <- DESeq(dds)
z<- results(res)
vsd <- vst(dds, blind=FALSE)
zz<-cbind(as.data.frame(z),assay(vsd))
dge<-as.data.frame(zz[order(zz$pvalue),])
write.table(dge,file="dge4_deseq.tsv",quote=F,sep="\t")
rnk<-as.data.frame( sign(dge$log2FoldChange) * (-log(dge$pvalue + 1E-307) ))
rownames(rnk)<-row.names(dge)
colnames(rnk)="Score"
rnk<-unique(merge(rnk,orth,by.x=0,by.y="Gene.stable.ID"))
rnk$Row.names=NULL
rnk<-aggregate(. ~ Human.gene.name,rnk,sum)
rnk<-rnk[-which(rnk$Human.gene.name==""),]
rownames(rnk)<-rnk$Human.gene.name
rnk$Human.gene.name=NULL
write.table(rnk,file="dge4.rnk",sep='\t',quote=F)
rnk4<-rnk
dge4<-dge

#some plots
pdf("dge4_plots.pdf")
sig<-subset(dge,padj<0.05)
SIG=nrow(sig)
DN=nrow(subset(sig,log2FoldChange<0))
UP=nrow(subset(sig,log2FoldChange>0))
HEADER=paste("HGMM vs HGWW:", SIG , "DGEs,", UP ,"upregulated,", DN, "downregulated")
plot(log2(dge$baseMean),dge$log2FoldChange,cex=0.6, xlab="log2 base mean", ylab="log2 fold change" ,pch=19,col="#838383")
points(log2(sig$baseMean),sig$log2FoldChange,cex=0.6,pch=19,col="red")
mtext(HEADER)
top<-head(sig,20)
#text(log2(top$baseMean)+1, top$log2FoldChange, labels = rownames(top),cex=0.7)
#volcano plot
plot(dge$log2FoldChange, -log2(dge$pvalue) ,cex=0.6, xlim=c(-4,6),xlab="log2 fold change", ylab="-log2 p-value" ,pch=19,col="#838383")
points(sig$log2FoldChange, -log2(sig$pvalue),cex=0.6,pch=19,col="red")
#text(top$log2FoldChange+0.5, -log2(top$pvalue), labels = rownames(top),cex=0.7)
mtext(HEADER)
# top N gene heatmap
colfunc <- colorRampPalette(c("blue", "white", "red"))
heatmap.2(  as.matrix(dge[1:100,c(7:ncol(dge))]), col=colfunc(25),scale="row",
 trace="none",margins = c(6,6), cexRow=.4, main="Top 100 genes")
dev.off()


