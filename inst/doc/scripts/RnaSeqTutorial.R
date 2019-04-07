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

##     # First uncompress your forward and reverse reads:

##     #!/bin/bash

##     mkdir ~/sortmerna

##    cp ~/share/Day01/data/sortmerna/* .

##     cd ~/sortmerna

##     mkdir ~/qa/sortmerna

##     mkdir ~/trimmomatic

##     cd ~/trimmomatic

##     mkdir ~/qa/trimmomatic

##     mkdir ~/qc_report


## ----child = "Chapters/05-Pseudo-Alignment.Rmd"--------------------------

##     cd

##     mkdir ~/kallisto


## ----child = "Chapters/12-R-Biological-QA.Rmd"---------------------------

## ----message=FALSE,warning=FALSE,results='hide'--------------------------
    library(RnaSeqTutorial)

## ------------------------------------------------------------------------
    library(tximport)
    files <- list.files("~/kallisto", 
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



