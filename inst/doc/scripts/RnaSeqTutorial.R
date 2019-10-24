## ----message=FALSE,warning=FALSE,results='hide',echo=FALSE---------------
options(digits=2)


## ----child = "Chapters/01-Convention.Rmd"--------------------------------

## ------------------------------------------------------------------------
a <- 1


## ------------------------------------------------------------------------
(a <- 1)



## ----child = "Chapters/02-Introduction.Rmd"------------------------------




## ----child = "Chapters/04-Preprocessing.Rmd"-----------------------------

##   # First, make an output directory called "qa" and a sub-directory of

##   # that called "raw"

##   mkdir -p ~/results/qa/raw

## 

##   # We can also make a directory to hold the raw data

##   mkdir ~/results/raw

## 

##   # Then link the files in our home folder. Always try to

##   # work on links or copies of the raw source data to keep them safe!

##   cd ~/results/raw

##   ln -s share/Day1/fastq/*.fq.gz .

## 

##   # Then you can run fastqc on the file(s) you have selected. FastQC can run

##   # a number of files in parallel (see the -t option)

##   fastqc -o ~/results/qa/raw -t 16 ~/results/raw/*.fq.gz

## 

##   # Finally, also run multiqc to summarize all data in one plot

##   cd ~/results/qa/raw

##   multiqc .


##     # THIS IS AN EXAMPLE OF THE COMMANDS! GO ON READING, do NOT execute them!

## 
##     # First uncompress your forward and reverse reads:

##     gunzip read_1.fq.gz read_2.fq.gz

## 
##     # Merge forward and reverse reads into a single file (this is a requirement

##     # of the sortmerna software)

##     merge-paired-reads.sh read_1.fq read_2.fq read-interleaved.fq

## 
##     # Run sortmerna

##     sortmerna --ref $SORTMERNA_DB --reads read-interleaved.fq --aligned \

##     read-rRNA-hits --other read-sortmerna --log -a 16 -v --paired_in --fastx

## 
##     # Split the file again in two separate files, one for forward, one for reverse

##     # orientation

##     unmerge-paired-reads.sh read-sortmerna.fq read-sortmerna_1.fq read-sortmerna_2.fq

## 
##     # Delete the interleaved files and compress the result files to save space

##     rm read-interleaved.fq read-sortmerna.fq

##     gzip read-sortmerna_1.fq read-sortmerna_2.fq


##     #!/bin/bash

## 

##     READ_FW="$1"

##     READ_RV="$2"

## 

##     FILEBASE=$(basename "${READ_FW/_1.fq.gz/}")

## 

##     echo "Uncompressing FASTQ data of $FILEBASE"

##     gunzip -c "$READ_FW" > ${FILEBASE}_1.fq

##     gunzip -c "$READ_RV" > ${FILEBASE}_2.fq

## 

##     READ_FW="${READ_FW%.gz}"

##     READ_RV="${READ_RV%.gz}"

## 

##     echo "Merging pairs of $FILEBASE"

##     merge-paired-reads.sh "$READ_FW" "$READ_RV" "${FILEBASE}_interleaved.fq"

## 

##     echo "Running SortMeRNA for $FILEBASE"

##     sortmerna --ref $SORTMERNA_DB --reads "${FILEBASE}_interleaved.fq" --aligned \

##     "${FILEBASE}-rRNA-hits" --other  "${FILEBASE}-sortmerna" --log -a 16 \

##     -v --paired_in --fastx

## 

##     echo "Unmerging SortMeRNA filtered pairs for $FILEBASE"

##     unmerge-paired-reads.sh "${FILEBASE}-sortmerna.fq" \

##     "${FILEBASE}-sortmerna_1.fq" "${FILEBASE}-sortmerna_2.fq"

## 

##     echo "Doing cleanup for $FILEBASE"

##     gzip "$READ_FW" "$READ_RV" "${FILEBASE}-sortmerna_1.fq" \

##     "${FILEBASE}-sortmerna_2.fq" "${FILEBASE}-rRNA-hits.fq"

##     rm "${FILEBASE}_interleaved.fq" "${FILEBASE}-sortmerna.fq"


##     mkdir ~/results/sortmerna

##     cd ~/results/sortmerna

##     find ~/results/raw -name "*.fq.gz"  | sort | head -n 4 | while read READ_FW

##     do

##       read READ_RV

##       bash ~/results/runSortMeRNA.sh $READ_FW $READ_RV

##     done


##    ln -s ~/share/Day1/sortmerna/* .


##     cd ~/results/sortmerna

##     multiqc .


##     mkdir ~/results/qa/sortmerna

##     fastqc -o ~/results/qa/sortmerna -t 16 ~/results/sortmerna/*sortmerna*.fq.gz


##     mkdir ~/results/trimmomatic

##     cd ~/results/trimmomatic

##     find ~/results/sortmerna -name "*sortmerna_[12].fq.gz" | sort | head -n 4 | while read FW_READ

##     do

##       read RV_READ

##       FILEBASE=$(basename "${FW_READ%_1.fq.gz}")

##       trimmomatic PE -threads 16 -phred64 "$FW_READ" "$RV_READ" \

##       "$FILEBASE-trimmomatic_1.fq.gz" "$FILEBASE-trimmomatic-unpaired_1.fq.gz" \

##       "$FILEBASE-trimmomatic_2.fq.gz" "$FILEBASE-trimmomatic-unpaired_2.fq.gz"  \

##       ILLUMINACLIP:"/usr/share/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa":2:30:10 \

##       SLIDINGWINDOW:5:20 MINLEN:50 2> "$FILEBASE-timmomatic.log"

##     done


##     cd ~/results/trimmomatic

##     ln -s ~/share/Day1/trimmomatic/* .

##     multiqc .


##     mkdir ~/results/qa/trimmomatic

##     fastqc -o ~/results/qa/trimmomatic -t 16 ~/results/trimmomatic/*trimmomatic_[12].fq.gz


##     mkdir ~/results/qc_report

##     ln -s ~/results/qa/raw/*zip ~/results/qc_report

##     ln -s ~/results/qa/sortmerna/*zip ~/results/qc_report

##     ln -s ~/results/qa/trimmomatic/*zip ~/results/qc_report

##     ln -s ~/results/sortmerna/*.log ~/results/qc_report

##     ln -s ~/results/trimmomatic/*.log ~/results/qc_report

##     multiqc -o ~/results/qc_report ~/results/qc_report



## ----child = "Chapters/05-Pseudo-Alignment.Rmd"--------------------------

##     cd ~/results

##     kallisto index -i Potra01-mRNA.idx \

##     ~/share/Day1/reference/fasta/Potra01-mRNA.fa.gz


##     mkdir -p ~/results/kallisto

##     cd ~/results/kallisto

##     find ../trimmomatic -name "*trimmomatic_[12].fq.gz" | sort | head -n 4 | while read FW_READ

##     do

##       read RV_READ

##       FILEBASE=$(basename "${FW_READ%_1.fq.gz}")

##       kallisto quant -i ~/results/Potra01-mRNA.idx -b 100 \

##       -o . -t 16 "$FW_READ" "$RV_READ" | tee $FILEBASE.log

##       # Kallisto doesn't let us specify an output filename so we rename all

##       # output files

##       mv "abundance.tsv" $FILEBASE"_abundance.tsv"

##       mv "abundance.h5" $FILEBASE"_abundance.h5"

##       mv "run_info.json" $FILEBASE"_run_info.json"

##     done


##     ~/results/kallisto

##     multiqc.


##     mkdir ~/results/qc_report

##     multiqc -o ~/results/qc_report ~/results



## ----child = "Chapters/12-R-Biological-QA.Rmd"---------------------------

## ----message=FALSE,warning=FALSE,results='hide'--------------------------
    library(RnaSeqTutorial)


## ------------------------------------------------------------------------
    library(tximport)
    files <- list.files("~/results/kallisto", 
                    pattern = ".*_abundance.tsv",
                    full.names = TRUE)
    tx <- suppressMessages(tximport(files = files,
                                     type = "kallisto", 
                                     txOut = TRUE))


## ------------------------------------------------------------------------
    tx.counts <- round(tx$counts)
    colnames(tx.counts) <- sub("_.*","",basename(files))


## ------------------------------------------------------------------------
tx2gene <- data.frame(
    TX=rownames(tx.counts),
    GENEID=sub("\\.\\d+$","",rownames(tx.counts)))


## ------------------------------------------------------------------------
    count.table <- round(summarizeToGene(tx,tx2gene)$counts)
    colnames(count.table) <- sub("_.*","",basename(files))


## ------------------------------------------------------------------------
sel <- rowSums(count.table) == 0
sprintf("%s percent of %s genes are not expressed",round(sum(sel) * 100/ nrow(count.table),digits=1),nrow(count.table))


## ------------------------------------------------------------------------
sel <- rowSums(tx.counts) == 0
sprintf("%s percent of %s transcripts are not expressed",round(sum(sel) * 100/ nrow(tx.counts),digits=1),nrow(tx.counts))


## ------------------------------------------------------------------------
library(RColorBrewer)
pal <- brewer.pal(8,"Dark2")
mar <- par("mar")
plot(density(log10(rowMeans(count.table))),col=pal[1],
     main="gene mean raw counts distribution",
     xlab="mean raw counts (log10)")


## ------------------------------------------------------------------------
plot.multidensity(log10(count.table),col=rep(pal,each=3),
                  legend.x="topright",legend.cex=0.5,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")


## ------------------------------------------------------------------------
plot.multidensity(log10(tx.counts),col=rep(pal,each=3),
                  legend.x="topright",legend.cex=0.5,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")


## ------------------------------------------------------------------------
samples <- read.csv(file.path(extdata(),"sex-samples.csv"))
sex <- samples$sex[match(colnames(count.table),samples$sample)]
date <- factor(samples$date[match(colnames(count.table),samples$sample)])


## ------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(sex=sex, date=date),
  design = ~ date + sex)


## ------------------------------------------------------------------------
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
sizes
boxplot(sizes,main="relative library sizes",ylab="scaling factor")


## ------------------------------------------------------------------------
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(count.table)


## ------------------------------------------------------------------------
vst <- vst - min(vst)


## ------------------------------------------------------------------------
meanSdPlot(assay(vsd)[rowSums(count.table)>0,])


## ------------------------------------------------------------------------
meanSdPlot(log2(as.matrix(count.table[rowSums(count.table)>0,])))


## ------------------------------------------------------------------------
meanSdPlot(log2(counts(dds,normalized=TRUE)[rowSums(count.table)>0,]))


## ------------------------------------------------------------------------
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)
smpls <- conditions


## ------------------------------------------------------------------------
sex.cols<-c("pink","lightblue")
sex.names<-c(F="Female",M="Male")


## ------------------------------------------------------------------------
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=sex.cols[as.integer(factor(sex))],
              pch=19)
legend("bottomright",pch=c(NA,15,15),col=c(NA,sex.cols[1:2]),
       legend=c("Color:",sex.names[levels(factor(sex))]))
par(mar=mar)


## ------------------------------------------------------------------------
dat <- date
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(factor(dat))],
              pch=19)
legend("topleft",pch=rep(c(19,23),each=10),col=rep(pal,2),legend=levels(factor(dat)),bty="n")
par(mar=mar)


## ------------------------------------------------------------------------
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=sex.cols[as.integer(factor(sex))],
              pch=c(19,17)[as.integer(factor(dat))])
legend("bottomright",pch=c(NA,15,15),col=c(NA,sex.cols[1:2]),
       legend=c("Color:",sex.names[levels(factor(sex))]))
legend("topleft",pch=c(NA,21,24),col=c(NA,1,1),
       legend=c("Symbol:",levels(factor(dat))),cex=0.85)
par(mar=mar)


## ------------------------------------------------------------------------
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=sex.cols[as.integer(factor(sex))],
     pch=c(19,17)[as.integer(factor(dat))],
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("bottomleft",pch=c(NA,15,15),col=c(NA,sex.cols[1:2]),
       legend=c("Color:",sex.names[levels(factor(sex))]))
legend("topright",pch=c(NA,21,24),col=c(NA,1,1),
       legend=c("Symbol:",levels(factor(dat))),cex=0.85)
text(pc$x[,1],  
     pc$x[,2],
     labels=colnames(count.table),cex=.5,adj=-1)


## ------------------------------------------------------------------------
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=sex.cols[as.integer(factor(sex))],
     pch=c(19,17)[as.integer(factor(dat))],     
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("bottomleft",pch=c(NA,15,15),col=c(NA,sex.cols[1:2]),
       legend=c("Color:",sex.names[levels(factor(sex))]))
legend("topright",pch=c(NA,21,24),col=c(NA,1,1),
       legend=c("Symbol:",levels(factor(dat))),cex=0.85)



## ----child = "Commons/17-SessionInfo.Rmd"--------------------------------

## ----session info, echo=FALSE--------------------------------------------
sessionInfo()



## ----child = "Commons/14-Acknowledgments.Rmd"----------------------------




## ----child = "Commons/15-Footnotes.Rmd"----------------------------------




## ----child = "Commons/16-Images.Rmd"-------------------------------------



