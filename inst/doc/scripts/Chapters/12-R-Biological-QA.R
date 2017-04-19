## ----echo=FALSE,eval=FALSE,message=FALSE,warning=FALSE,results='hide'----
##     # TODO keep in mind a comparison of HTSeq with 11-R-Counts

## ----message=FALSE,warning=FALSE,results='hide'--------------------------
    library(RnaSeqTutorial)

## ------------------------------------------------------------------------
res <- mclapply(dir(file.path(extdata(),"htseq"),
                    pattern="^[2,3].*_STAR\\.txt",
                    full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=2)
names(res) <- gsub("_.*_STAR\\.txt","",dir(file.path(extdata(),"htseq"),
                                           pattern="^[2,3].*_STAR\\.txt"))

## ------------------------------------------------------------------------
addInfo <- c("__no_feature","__ambiguous",
             "__too_low_aQual","__not_aligned",
             "__alignment_not_unique")

## ------------------------------------------------------------------------
sel <- match(addInfo,res[[1]][,1])
count.table <- do.call(cbind,lapply(res,"[",2))[-sel,]
colnames(count.table) <- names(res)
rownames(count.table) <- res[[1]][,1][-sel]

## ------------------------------------------------------------------------
count.stats <- do.call(cbind,lapply(res,"[",2))[sel,]
colnames(count.stats) <- names(res)
rownames(count.stats) <- sub("__","",addInfo)
count.stats <- rbind(count.stats,aligned=colSums(count.table))
count.stats <- count.stats[rowSums(count.stats) > 0,]

## ------------------------------------------------------------------------
apply(count.stats,2,function(co){round(co*100/sum(co))})

## ------------------------------------------------------------------------
pal=brewer.pal(6,"Dark2")[1:nrow(count.stats)]
mar <- par("mar")
par(mar=c(7.1,5.1,4.1,2.1))
barplot(as.matrix(count.stats),col=pal,beside=TRUE,las=2,main="read proportion",
        ylim=range(count.stats)+c(0,1e+7))
legend("top",fill=pal,legend=rownames(count.stats),bty="n",cex=0.8,horiz=TRUE)
par(mar=mar)

## ------------------------------------------------------------------------
sel <- rowSums(count.table) == 0
sprintf("%s percent",round(sum(sel) * 100/ nrow(count.table),digits=1))
sprintf("of %s genes are not expressed",nrow(count.table))

## ------------------------------------------------------------------------
plot(density(log10(rowMeans(count.table))),col=pal[1],
     main="mean raw counts distribution",
     xlab="mean raw counts (log10)")

## ------------------------------------------------------------------------
pal=brewer.pal(8,"Dark2")
plot.multidensity(log10(count.table),col=rep(pal,each=3),
                  legend.x="topright",legend.cex=0.5,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

## ------------------------------------------------------------------------
samples <- read.csv(file.path(extdata(),"sex-samples.csv"))
conditions <- names(res)
sex <- samples$sex[match(names(res),samples$sample)]
date <- factor(samples$date[match(names(res),samples$sample)])

## ------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=conditions,
                       sex=sex,
                       date=date),
  design = ~ condition)

## ------------------------------------------------------------------------
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- names(res)
sizes
boxplot(sizes,main="relative library sizes",ylab="scaling factor")

## ------------------------------------------------------------------------
colData(dds)$condition <- factor(colData(dds)$condition,
                                 levels=unique(conditions))
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
     labels=conditions,cex=.5,adj=-1)

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

## ------------------------------------------------------------------------
sel <- order(apply(vst,1,sd),decreasing=TRUE)[1:1000]

## ------------------------------------------------------------------------
heatmap.2(vst[sel,],labRow = NA,trace = "none",cexCol = 0.6 )

## ------------------------------------------------------------------------
heatmap.2(vst[sel,],labRow = NA, labCol = sex,trace = "none",cexCol = 0.6 )

## ------------------------------------------------------------------------
heatmap.2(vst[sel,],labRow = NA, labCol = date,trace = "none",cexCol = 0.6 )

