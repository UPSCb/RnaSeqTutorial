## ----message=FALSE,warning=FALSE,results='hide',echo=FALSE---------------
    options(digits=2)
    library(RnaSeqTutorial)
    bamfiles <- dir(file.path(extdata(),"BAM"),
                   pattern="*.bam$",full.names=TRUE)

## ------------------------------------------------------------------------
    bwtFiles <- dir(path=file.path(extdata(),"bowtie"),
                   pattern="*.bwt$",full.names=TRUE)

## ----eval=FALSE----------------------------------------------------------
##     read.delim(file=file.path(extdata(),"bowtie","SRR074430.bwt"),
##               header=FALSE,nrows=10)

## ----eval=FALSE----------------------------------------------------------
##     alignedRead <- readAligned(dirPath=file.path(extdata(),"bowtie"),
##                        pattern="*.bwt$",type="Bowtie")

## ----eval=FALSE----------------------------------------------------------
##     alignedRead2L <- readAligned(dirPath=file.path(extdata(),"bowtie"),
##                                 pattern="*.bwt$",type="Bowtie",
##                        filter=chromosomeFilter("2L"))

## ------------------------------------------------------------------------
    fl <- system.file("extdata", "ex1.sam", package="Rsamtools")
    strsplit(readLines(fl, 1), "\t")[[1]]

## ------------------------------------------------------------------------
    alnFile <- system.file("extdata", "ex1.bam", package="Rsamtools")
    aln <- readGAlignments(alnFile)

## ----eval=FALSE----------------------------------------------------------
##     head(aln, 3)

## ------------------------------------------------------------------------
    table(strand(aln))
    table(width(aln))
    head(sort(table(cigar(aln)), decreasing=TRUE))

## ---- eval=FALSE---------------------------------------------------------
##     param <- ScanBamParam(what="seq")
##     seqs <- scanBam(bamfiles[[1]], param=param)
##     readGC <- gcFunction(seqs[[1]][["seq"]])
##     hist(readGC)

## ----eval=FALSE----------------------------------------------------------
##     scanBamHeader(bamfiles[1])

## ---- eval=FALSE---------------------------------------------------------
##     param <- ScanBamParam(tag="NH")
##     nhs <- scanBam(bamfiles[[1]], param=param)[[1]]$tag$NH

## ----eval=FALSE----------------------------------------------------------
##     param <- ScanBamParam(what="cigar")
##     cigars <- scanBam(bamfiles[[1]], param=param)[[1]]$cigar
##     cigar.matrix <- cigarOpTable(cigars)
##     intron.size <- cigar.matrix[,"N"]
##     intron.size[intron.size>0]
##     plot(density(intron.size[intron.size>0]))
##     hist(log10(intron.size[intron.size>0]),xlab="intron size (log10 bp)")

## ----eval=FALSE----------------------------------------------------------
##     library(leeBamViews)
##     bpaths = dir(system.file("bam", package="leeBamViews"), full=TRUE, patt="bam$")
##     gt<-do.call(rbind,strsplit(basename(bpaths),"_"))[,1]
##     geno<-substr(gt,1,nchar(gt)-1)
##     lane<-substr(gt,nchar(gt),nchar(gt))
##     pd = DataFrame(geno=geno, lane=lane, row.names=paste(geno,lane,sep="."))
##     bs1 = BamViews(bamPaths=bpaths, bamSamples=pd,
##                    bamExperiment=list(annotation="org.Sc.sgd.db"))
##     bamPaths(bs1)
##     bamSamples(bs1)

## ----eval=FALSE----------------------------------------------------------
##     sel <- GRanges(seqnames = "Scchr13",
##                    IRanges(start = 861250, end = 863000),
##                    strand="+")
##     covex = RleList(lapply(bamPaths(bs1),
##                            function(x){
##                              coverage(readGAlignments(x))[sel][["Scchr13"]]}))

