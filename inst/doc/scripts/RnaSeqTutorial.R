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

##    cp ~/share/Day02/data/sortmerna/* .

##     cd ~/sortmerna

##     mkdir ~/qa/sortmerna

##     mkdir ~/trimmomatic

##     cd ~/trimmomatic

##     mkdir ~/qa/trimmomatic

##     mkdir ~/qc_report


## ----child = "Chapters/05-Pseudo-Alignment.Rmd"--------------------------

##     cd

##     mkdir ~/kallisto


## ----child = "Chapters/14-R-Biological-QA-FG.Rmd"------------------------

## ----message=FALSE,warning=FALSE,results='hide'--------------------------
    library(RnaSeqTutorial)

## ------------------------------------------------------------------------
    counts <- readAbundance("~/share/sex/kallisto")

## ------------------------------------------------------------------------
    nonExpressed(counts)

## ------------------------------------------------------------------------
rawDataMeanPlot(counts)

## ------------------------------------------------------------------------
rawDataSamplePlot(counts)

## ------------------------------------------------------------------------
samples <- read.csv("~/share/sex/doc/samples.csv")
dds <- createDESeqDataSet(counts,samples)

## ------------------------------------------------------------------------
reportSizeFactors(dds)

## ------------------------------------------------------------------------
vst <- transform(dds)

## ------------------------------------------------------------------------
vst <- vst - min(vst)

## ------------------------------------------------------------------------
validateVST(vst)

## ------------------------------------------------------------------------
plotUnTransformed(dds)

## ------------------------------------------------------------------------
plotPca(vst,samples)


## ----child = "Chapters/15-R-Differential-Expression-FG.Rmd"--------------

## ----echo=FALSE,eval=FALSE,message=FALSE,warning=FALSE,results='hide'----
##     # TODO keep in mind to save (as an rda in data) a copy of the count table

## ----message=FALSE,warning=FALSE,results='hide'--------------------------
    library(RnaSeqTutorial)
    counts <- readAbundance("~/share/sex/kallisto")
    samples <- read.csv("~/share/sex/doc/samples.csv")
    dds <- createDESeqDataSet(counts,samples)

## ------------------------------------------------------------------------
# the object
dds

# the metadata
colData(dds)

# the design
design(dds)

## ------------------------------------------------------------------------
dds <- DESeq(dds)

## ------------------------------------------------------------------------
resultsNames(dds)

## ------------------------------------------------------------------------
alpha=0.01
log2FC=0.5
res <- results(dds,name = "sex_M_vs_F")

## ------------------------------------------------------------------------
res

## ------------------------------------------------------------------------
res[abs(res$log2FoldChange) >= log2FC &
        ! is.na(res$padj) &
        res$padj <= alpha,]

## ------------------------------------------------------------------------
plotMA(as(res,"DataFrame"),
       alpha=alpha,
       log2FC=log2FC)

## ------------------------------------------------------------------------
volcanoPlot(res,alpha=alpha,log2FC=log2FC)

## ------------------------------------------------------------------------
# Potri.019G047300, grouping by sex
dotplot(counts(dds,normalized=TRUE)["Potri.019G047300",],
        groups=samples$sex,col=c("pink","lightblue"),pch=19,
        xlab="library size corrected counts")

# Potri.019G047300, grouping by date
dotplot(counts(dds,normalized=TRUE)["Potri.019G047300",],
        groups=samples$date,col=c("darkgreen","brown"),pch=19,
                xlab="library size corrected counts")


## ------------------------------------------------------------------------
# Potri.014G155300, grouping by sex
dotplot(counts(dds,normalized=TRUE)["Potri.014G155300",],
        groups=samples$sex,col=c("pink","lightblue"),pch=19,
        xlab="library size corrected counts")

# Potri.014G155300, grouping by date
dotplot(counts(dds,normalized=TRUE)["Potri.014G155300",],
        groups=samples$date,col=c("darkgreen","brown"),pch=19,
                xlab="library size corrected counts")


## ----child = "Chapters/16-R-DE-practical.Rmd"----------------------------



## ----child = "Commons/17-SessionInfo.Rmd"--------------------------------

## ----session info, echo=FALSE--------------------------------------------
sessionInfo()


## ----child = "Commons/14-Acknowledgments.Rmd"----------------------------



## ----child = "Commons/18-Appendix.Rmd"-----------------------------------

## ------------------------------------------------------------------------
showMethods(readAbundance,includeDefs = TRUE)

## ------------------------------------------------------------------------
showMethods(nonExpressed,includeDefs = TRUE)


## ----child = "Commons/15-Footnotes.Rmd"----------------------------------



## ----child = "Commons/16-Images.Rmd"-------------------------------------



