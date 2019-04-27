#!/bin/bash

run_gsea(){
GSEAJAR=*jar
RNK=$1
GMT=$2
echo $GSEAJAR $RNK $GMT
java -Xmx4096m -cp $GSEAJAR xtools.gsea.GseaPreranked  \
-gmx $GMT -collapse false -mode Max_probe \
-norm meandiv -nperm 1000 -rnk $RNK -scoring_scheme classic \
-rpt_label ${RNK}.${GMT} -include_only_symbols true -make_sets true \
-plot_top_x 20 -rnd_seed timestamp -set_max 10000 -set_min 10 -zip_report false \
-out . -gui false
}
export -f run_gsea



sed -i 's/Score/Accession_GeneID\tScore/' *rnk
for RNK in *rnk ; do
  sed 1d $RNK | cut -d '_' -f2- | awk '{OFS="\t"} {print $1,$2,$2*$2}' \
  | sort -k3gr | awk '{OFS="\t"} !arr[$1]++ {print $1,$2}' > tmp
  sed -e '1i\GeneID\tScore' tmp > $RNK
done


parallel -j6 run_gsea ::: *.rnk ::: *.gmt

echo 'GeneSetName GeneSetSize ES NES p-val FDR FWER' | sed 's/ /\t/g' > header.txt

for GSEADIR in `ls | grep GseaPreranked | grep -v xls$` ; do
  awk '{FS="\t"} {OFS="\t"} $8<0.05 {print $1,$4,$5,$6,$7,$8,$9} ' $GSEADIR/gsea_report_for_na_*xls \
  | cat header.txt - > $GSEADIR.xls
done

cat <<'EOF' >> plotGSEA.R
path = "."
infiles <- dir(pattern='\\.xls$')
plot.GSEA <- function(file){
  x<-fread(file)
  rownames(x)<-mytest$GeneSetName
  x <-read.table(file,header=T,row.names=1)
  y<-head(x[order(x$ES),],n=20L)
  y<-y[which(y$ES<0),]
  z<-head(x[order(-x$ES),],n=20L)
  z<-z[which(z$ES>0),]
  df <- rbind(y,z)
  df<-df[order(df$ES),]
  barplot(df$ES,main=file,xlab="GSEA enrichment score",horiz=TRUE,names.arg=row.names(df))
}
pdf(file="gsea_results.pdf",width=15,height=10)
par(las=2) ; par(mar=c(10,50,1,1))
lapply(infiles,plot.GSEA)
dev.off()
lapply(infiles,plot.GSEA)
EOF
Rscript plotGSEA.R

