## ----echo=FALSE,eval=FALSE,message=FALSE,warning=FALSE,results='hide'----
##     # TODO keep in mind to save (as an rda in data) a copy of the count table

## ----message=FALSE,warning=FALSE,results='hide'--------------------------
    library(RnaSeqTutorial)
    samples <- read.csv(file.path(extdata(),"sex-samples.csv"))
    res <- mclapply(dir(file.path(extdata(),"htseq"),
                    pattern="^[2,3].*_STAR\\.txt",
                    full.names=TRUE),function(fil){
      read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
    },mc.cores=2)
    names(res) <- gsub("_.*_STAR\\.txt","",dir(file.path(extdata(),"htseq"),
                                           pattern="^[2,3].*_STAR\\.txt"))
    addInfo <- c("__no_feature","__ambiguous",
             "__too_low_aQual","__not_aligned",
             "__alignment_not_unique")
    sel <- match(addInfo,res[[1]][,1])
    count.table <- do.call(cbind,lapply(res,"[",2))[-sel,]
    colnames(count.table) <- names(res)
    rownames(count.table) <- res[[1]][,1][-sel]
    conditions <- names(res)
    sex <- samples$sex[match(names(res),samples$sample)]
    date <- factor(samples$date[match(names(res),samples$sample)])
    pal=brewer.pal(8,"Dark2")

## ------------------------------------------------------------------------
yearSexDesign <- samples[grep("[2,3].*",samples$sample),c("sex","date")]

## ------------------------------------------------------------------------
cdsFull <- newCountDataSet(count.table,yearSexDesign)
cdsFull = estimateSizeFactors( cdsFull )
cdsFull = estimateDispersions( cdsFull )

## ------------------------------------------------------------------------
plotDispLSD(cdsFull)

## ------------------------------------------------------------------------
fit1 = suppressWarnings(fitNbinomGLMs( cdsFull, count ~ date + sex ))
fit0 = suppressWarnings(fitNbinomGLMs( cdsFull, count ~ date))

## ------------------------------------------------------------------------
sel <- !is.na(fit1$converged) & !is.na(fit0$converged) & fit1$converged & fit0$converged
fit1 <- fit1[sel,]
fit0 <- fit0[sel,]

## ------------------------------------------------------------------------
pvalsGLM = nbinomGLMTest( fit1, fit0 )
padjGLM = p.adjust( pvalsGLM, method="BH" )

## ------------------------------------------------------------------------
boxplot(rowSums(count.table[rownames(fit1[fit1$sexM > 20,]),]>0),
        ylab="number of sample with expression > 0",
        main="distribution of the # of samples with\nexpression > 0 for genes with extreme high log2FC"
)
boxplot(rowSums(count.table[rownames(fit1[fit1$sexM < -20,]),]>0),
        ylab="number of sample with expression > 0",
        main="distribution of the # of samples with\nexpression > 0 for genes with extreme low log2FC"
)

## ------------------------------------------------------------------------
plot(density(padjGLM[fit1$sexM > 20]),
     main="Adj. p-value for genes with extreme log2FC",
     xlab="Adj. p-value")
plot(density(padjGLM[fit1$sexM > -20]),
     main="Adj. p-value for genes with extreme log2FC",
     xlab="Adj. p-value")

## ------------------------------------------------------------------------
sel <- fit1$sexM > -20 & fit1$sexM < 20
fit1 <- fit1[sel,]
padjGLM <- padjGLM[sel]
pvalsGLM <- pvalsGLM[sel]

## ------------------------------------------------------------------------
padjGLM[padjGLM==0] <- min(padjGLM[padjGLM!=0])/10
pvalsGLM[pvalsGLM==0] <- min(pvalsGLM[pvalsGLM!=0])/10

## ------------------------------------------------------------------------
heatscatter(fit1$sexM,-log10(pvalsGLM),
            main="Male vs. Female Differential Expression",cor=FALSE,
            xlab="Log2 Fold Change", ylab="- log(10) p-value")
legend("topleft",bty="n","1% p-value cutoff",lty=2,col="gray")
points(fit1[pvalsGLM<.01,"sexM"],
       -log10(pvalsGLM[pvalsGLM<.01]),
       col="lightblue",pch=19)
pos <- c(order(pvalsGLM)[1:4],
         order(fit1$sexM,decreasing=TRUE)[1:4],
         order(fit1$sexM)[1:4])
points(fit1[pos,"sexM"],-log10(pvalsGLM[pos]),col="red")
abline(h=2,lty=2,col="gray")

## ------------------------------------------------------------------------
heatscatter(fit1$sexM,-log10(padjGLM),
            main="Male vs. Female Differential Expression",cor=FALSE,
            xlab="Log2 Fold Change", ylab="- log(10) adj. p-value")
legend("topleft",bty="n","1% FDR cutoff",lty=2,col="gray")
points(fit1[padjGLM<.01,"sexM"],-log10(padjGLM[padjGLM<.01]),col="lightblue",pch=19)
pos <- c(order(padjGLM)[1:4],order(fit1$sexM,decreasing=TRUE)[1:4],order(fit1$sexM)[1:4])
points(fit1[pos,"sexM"],-log10(padjGLM[pos]),col="red")
abline(h=2,lty=2,col="gray")

## ------------------------------------------------------------------------
UmAspD <- rownames(fit1[padjGLM<.01,])

## ------------------------------------------------------------------------
UmAspD

## ------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
    countData = count.table,
    colData = data.frame(condition=conditions,
                       sex=sex,
                       date=date),
    design = ~date+sex)
dds <- DESeq(dds)
plotDispEsts(dds)
res <- results(dds,contrast=c("sex","M","F"))

## ------------------------------------------------------------------------
alpha=0.01
plotMA(res,alpha=alpha)
volcanoPlot(res,alpha=alpha)

## ------------------------------------------------------------------------
hist(res$padj,breaks=seq(0,1,.01))
sum(res$padj<alpha,na.rm=TRUE)

## ------------------------------------------------------------------------
UmAspD2 <- rownames(res[order(res$padj,na.last=TRUE),][1:sum(res$padj<alpha,na.rm=TRUE),])
UmAspD2

## ------------------------------------------------------------------------
as.data.frame(res)["Potra000503g03273.0",]
data.frame(
  counts=counts(dds,normalized=TRUE)["Potra000503g03273.0",],
  sex,
  date,
  conditions)
mcols(dds)[names(rowRanges(dds))=="Potra000503g03273.0",]

## ------------------------------------------------------------------------
sex[conditions=="226.1"] <- "M"
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=conditions,
                       sex=sex,
                       date=date),
  design = ~ date+sex)
dds <- DESeq(dds)
res.swap <- results(dds,contrast=c("sex","M","F"))
plotMA(res.swap,alpha=alpha)
volcanoPlot(res.swap,alpha=alpha)
hist(res.swap$padj,breaks=seq(0,1,.01))
sum(res.swap$padj<alpha,na.rm=TRUE)
UmAspD2.swap <- rownames(res.swap[order(res.swap$padj,na.last=TRUE),][1:sum(res.swap$padj<alpha,na.rm=TRUE),])
UmAspD2.swap

## ------------------------------------------------------------------------
sel <- conditions != "226.1"
dds <- DESeqDataSetFromMatrix(
  countData = count.table[,sel],
  colData = data.frame(condition=conditions[sel],
                       sex=sex[sel],
                       date=date[sel]),
  design = ~ date+sex)
dds <- DESeq(dds)
plotDispEsts(dds)
res.sel <- results(dds,contrast=c("sex","M","F"))
plotMA(res.sel,alpha=alpha)
volcanoPlot(res.sel,alpha=alpha)
hist(res.sel$padj,breaks=seq(0,1,.01))
sum(res.sel$padj<alpha,na.rm=TRUE)
UmAspD2.sel <- rownames(res.sel[order(res.sel$padj,na.last=TRUE),][1:sum(res.sel$padj<alpha,na.rm=TRUE),])
UmAspD2.sel

## ------------------------------------------------------------------------
plot.new()
grid.draw(venn.diagram(list(
  UmAsp=UmAspD2,
  UmAsp.swap=UmAspD2.swap,
  UmAsp.sel=UmAspD2.sel,
  UmAspD = UmAspD),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("UmAsp - DESeq2",
                                        "UmAsp - corrected",
                                        "UmAsp - removed",
                                        "UmAsp - DESeq")))

## ------------------------------------------------------------------------
sort(intersect(UmAspD2,UmAspD))
sort(intersect(UmAspD2.swap,UmAspD))
sort(intersect(UmAspD2.sel,UmAspD))

## ------------------------------------------------------------------------
## mis-classified
# sex <- samples$sex[match(colnames(count.table),samples$sample)]
# date <- factor(samples$date[match(colnames(count.table),samples$sample)])
# ddsM <- DESeqDataSetFromMatrix(
#   countData = count.table,
#   colData = data.frame(condition=conditions,
#                        sex=sex,
#                        date=date),
#   design = ~ date+sex)
# ddsM <- DESeq(ddsM)
# resM <- results(ddsM,contrast=c("sex","F","M"))
# 
# counts(ddsM,normalized=TRUE)["Potri.019G047300.0",]
# resM["Potri.019G047300.0",]
# mcols(ddsM)[names(rowData(ddsM)) == "Potri.019G047300.0",]
# 
# ## reclassified
# sexR <- sex
# sexR[5] <- "M"
# ddsR <- DESeqDataSetFromMatrix(
#   countData = count.table,
#   colData = data.frame(condition=conditions,
#                        sex=sexR,
#                        date=date),
#   design = ~ date+sex)
# ddsR <- DESeq(ddsR)
# resR <- results(ddsR,contrast=c("sex","F","M"))
# 
# counts(ddsR,normalized=TRUE)["Potri.019G047300.0",]
# resR["Potri.019G047300.0",]
# mcols(ddsR)[names(rowData(ddsR)) == "Potri.019G047300.0",]
# 
# ## samples details
# sex
# date
# 
# ## the model is not be affected by the reclassification
# lapply(split(sex,date),table)
# lapply(split(sexR,date),table)
# 

