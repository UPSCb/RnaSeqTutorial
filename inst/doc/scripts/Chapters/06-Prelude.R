## ----eval=FALSE----------------------------------------------------------
##    library(GenomicAlignments)
##     bamfiles <- dir(file.path(extdata(),"BAM"),
##                     pattern="*.bam$",full.names=TRUE)
##     gAlns <- lapply(bamfiles,readGAlignments)

## ----eval=FALSE----------------------------------------------------------
##     library(parallel)
##     gAlns <- mclapply(bamfiles,readGAlignments)

## ----eval=FALSE----------------------------------------------------------
##     library(Rsamtools)
##     bamFileList <- BamFileList(bamfiles,yieldSize=10^5)

## ----eval=FALSE----------------------------------------------------------
##     ## DO NOT RUN ME!
##     ## I AM JUST AN EXAMPLE
##     open(bamFile)
##     while(length(chunk <- readGAlignmentsFromBam(bamFile))){
##           message(length(chunk))
##     }
##     close(bamFile)

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

