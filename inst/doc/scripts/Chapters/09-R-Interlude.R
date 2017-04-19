## ----message=FALSE,warning=FALSE,results='hide'--------------------------
    library(RnaSeqTutorial)
    bamfiles <- dir(file.path(extdata(),"BAM"),
                   pattern="*.bam$",full.names=TRUE)
    bamFileList <- BamFileList(bamfiles,yieldSize=10^6)

## ----eval=FALSE----------------------------------------------------------
##     gAlns <- mclapply(bamFileList,function(bamFile){
##      open(bamFile)
##      gAln <- GAlignments()
##      while(length(chunk <- readGAlignments(bamFile))){
##        gAln <- c(gAln,chunk)
##      }
##      close(bamFile)
##      return(gAln)
##     })

## ------------------------------------------------------------------------
    data("Ptrichocarpa_v3.0_210_synthetic_transcripts")
    count.list <- mclapply(bamFileList,function(bamFile,annot){
    open(bamFile)
    counts <- vector(mode="integer",length=length(annot))
    while(length(chunk <- readGAlignments(bamFile))){
       counts <- counts + assays(summarizeOverlaps(annot,
                                                   chunk,
                                                   mode="Union",
                                                   ignore.strand=TRUE))$counts
    }
    close(bamFile)
    return(counts)
    },syntheticTrx[syntheticTrx$type == "exon"])

## ------------------------------------------------------------------------
    count.table <- do.call("cbind",count.list)
    head(count.table[rowSums(count.table)>0,])

