#' Functional Genomics High-Level functions
#'
#' High level functions for the Functional Genomics course computer lab
#'
#' @name Functional Genomics High Level
#' @rdname RnaSeqTutorial-functionalGenomics
#' @aliases readAbundance readAbundance,character-method nonExpressed
#' nonExpressed,matrix-method rawDataMeanPlot rawDataMeanPlot,matrix-method
#' rawDataSamplePlot rawDataSamplePlot,matrix-method createDESeqDataSet
#' createDESeqDataSet,matrix,data.frame-method reportSizeFactors
#' reportSizeFactors,DESeqDataSet-method transform transform,DESeqDataSet-method
#' plotUnTransformed plotUnTransformed,DESeqDataSet-method plotPca,
#' plotPca,matrix,data.frame-method
#' @param x a character value (readAbundance), a DESeqDataSet (reportSizeFactor,
#' transform) or a matrix (all other functions)
#' @param y a data.frame (createDESeqDataSet amd plotPca)
#' @param pattern a character pattern to match the abundance files (readAbundance)
#' @param type a character value to define the type of quantification (readAbundance)
#' @param design a design expression (createDESeqDataSet)
#' @return TODO FIXME
#' @examples
#' \dontrun{
#'   counts <- readAbundance(x="path-to-kallisto-output-dir")
#'   nonExpressed(counts)
#'   rawDataMeanPlot(counts)
#'   rawDataSamplePlot(counts)
#' }

# Generics
#' @exportMethod readAbundance
setGeneric(name="readAbundance",
           def=function(x,
                        pattern=character(0),
                        type=c("kallisto","salmon")){
               standardGeneric("readAbundance")
           })

#' @exportMethod nonExpressed
setGeneric(name="nonExpressed",
           def=function(x){
               standardGeneric("nonExpressed")
           })

#' @exportMethod rawDataMeanPlot
setGeneric(name="rawDataMeanPlot",
           def=function(x){
               standardGeneric("rawDataMeanPlot")
           })

#' @exportMethod rawDataSamplePlot
setGeneric(name="rawDataSamplePlot",
           def=function(x){
               standardGeneric("rawDataSamplePlot")
           })

#' @exportMethod createDESeqDataSet
setGeneric(name="createDESeqDataSet",
           def=function(x,y,design=~date+sex){
               standardGeneric("createDESeqDataSet")
           })

#' @exportMethod reportSizeFactors
setGeneric(name="reportSizeFactors",
           def=function(x){
               standardGeneric("reportSizeFactors")
           })


#' @exportMethod transform
setGeneric(name="transform",
           def=function(x){
               standardGeneric("transform")
           })

#' @exportMethod validateVST
setGeneric(name="validateVST",
           def=function(x){
               standardGeneric("validateVST")
           })

#' @exportMethod plotUnTransformed
setGeneric(name="plotUnTransformed",
           def=function(x){
               standardGeneric("plotUnTransformed")
           })

#' @exportMethod plotPca
setGeneric(name="plotPca",
           def=function(x,y,color=NULL,...){
               standardGeneric("plotPca")
           })

# Implementation
setMethod(f="readAbundance",
          signature="character",
          definition=function(x,
                              pattern=".*_abundance.tsv",
                              type=c("kallisto","salmon")){

              stopifnot(require(tximport))

              type <- match.arg(type)

              if(length(x)!=1){
                  stop("The 'x' argument needs to be an existing directory")
              }

              files <- list.files(x,
                                  pattern = pattern,
                                  recursive = TRUE,
                                  full.names = TRUE)

              if(length(files)==0){
                  stop("The directory provided as argument 'x' does not contain abundance files.")
              }

              # we read in the files
              tx <- suppressMessages(tximport(files = files,
                                              type = type,
                                              txOut = TRUE))

              # create the transcript to gene mapping
              tx2gene <- data.frame(
                  TX=rownames(tx$counts),
                  GENEID=sub("\\.\\d+$","",rownames(tx$counts)))

              # summarising gene expression
              count.table <- round(summarizeToGene(tx,tx2gene)$counts)
              gnames <- rownames(count.table)
              count.table <- apply(count.table,2,as.integer)
              rownames(count.table) <- gnames
              colnames(count.table) <- sub("_.*","",basename(files))

              return(count.table)
          })


setMethod(f="nonExpressed",
          signature="matrix",
          definition=function(x){

              sel <- rowSums(x) == 0
              message(sprintf("%s percent of %s genes are not expressed",
                      round(sum(sel) * 100/ nrow(x),digits=1),
                      nrow(x)))
})

setMethod(f="rawDataMeanPlot",
          signature="matrix",
          definition=function(x){

              require(RColorBrewer)
              pal <- brewer.pal(8,"Dark2")

              invisible(plot(density(log10(rowMeans(x))),col=pal[1],
                   main="gene mean raw counts distribution",
                   xlab="mean raw counts (log10)",lwd=2))

          })

setMethod(f="rawDataSamplePlot",
          signature="matrix",
          definition=function(x){

              require(RColorBrewer)
              pal <- brewer.pal(8,"Dark2")

              invisible(plot.multidensity(log10(x),col=rep(pal,each=3),
                                legend.x="topright",legend.cex=0.5,
                                main="sample raw counts distribution",
                                xlab="per gene raw counts (log10)"))
          })

setMethod(f="createDESeqDataSet",
          signature=c("matrix","data.frame"),
          definition=function(x,y,design= ~ date + sex){

              require(DESeq2)
              DESeqDataSetFromMatrix(
                  countData = x,
                  colData = y,
                  design = design)
          })

setMethod(f="reportSizeFactors",
          signature="DESeqDataSet",
          definition=function(x){
              require(DESeq2)
              sizes <- sizeFactors(estimateSizeFactors(x))
              boxplot(sizes,main="relative library sizes",ylab="scaling factor")
          })

setMethod(f="transform",
          signature="DESeqDataSet",
          definition=function(x){
              vst <- assay(varianceStabilizingTransformation(x, blind=TRUE))
              vst <- vst - min(vst)
          })

setMethod(f="validateVST",
          signature="matrix",
          definition=function(x){
              require(hexbin)
              require(vsn)
              meanSdPlot(x[rowSums(x)>0,])
          })

setMethod(f="plotUnTransformed",
          signature="DESeqDataSet",
          definition=function(x){
              sel <- rowSums(assay(dds))>0
              suppressWarnings(meanSdPlot(log2(counts(x,normalized=FALSE)[sel,])))
              dds <- estimateSizeFactors(x)
              suppressWarnings(meanSdPlot(log2(counts(dds,normalized=TRUE)[sel,])))
})

setMethod(f="plotPca",
          signature=c("matrix","data.frame"),
          definition=function(x,y,
                              color=NULL,...){

              pc <- prcomp(t(x))
              percent <- round(summary(pc)$importance[2,]*100)

              if(is.null(color)){
                  sex.cols<-c("pink","lightblue")
                  sex.names<-c(F="Female",M="Male")
                  col=sex.cols[as.integer(y$sex)]
                  symbol=c(19,17)[as.integer(y$date)]
              } else {
                  col=color
                  symbol=19
              }

              plot(pc$x[,1],
                   pc$x[,2],
                   xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
                   ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
                   col=col,
                   pch=symbol,
                   main="Principal Component Analysis",
                   sub="variance stabilized counts",
                   ...)

              if(is.null(color)){
                  legend("bottomleft",pch=c(NA,15,15),col=c(NA,sex.cols[1:2]),
                         legend=c("Color:",sex.names[levels(y$sex)]))

                  legend("topright",pch=c(NA,21,24),col=c(NA,1,1),
                         legend=c("Symbol:",sub("Y","",levels(y$dat))),cex=0.85)
              }

              text(pc$x[,1],
                   pc$x[,2],
                   labels=colnames(x),cex=.5,adj=-0.5)
          })
