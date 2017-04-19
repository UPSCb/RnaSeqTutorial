## ----message=FALSE,warning=FALSE,results='hide'--------------------------
    library(RnaSeqTutorial)
    bamfiles <- dir(file.path(extdata(),"BAM"),
                   pattern="*.bam$",full.names=TRUE)

## ------------------------------------------------------------------------
    bamfiles <- dir(file.path(extdata(), "BAM"), ".bam$", full=TRUE)
    names(bamfiles) <- sub("_.*", "", basename(bamfiles))

## ------------------------------------------------------------------------
    aln <- readGAlignments(bamfiles[1])
    strand(aln) <- "*"

## ------------------------------------------------------------------------
    param <- ScanBamParam(tag="NH")
    nhs <- scanBam(bamfiles[[1]], param=param)[[1]]$tag$NH
    aln <- aln[nhs==1,]

## ----eval=FALSE----------------------------------------------------------
##     load("~/Ptrichocarpa_v3.0_210_synthetic_transcripts.rda")
##     annot <- syntheticTrx[syntheticTrx$type=="mRNA",]
##     counts1 <- summarizeOverlaps(annot, aln, mode="Union")
##     counts2 <- summarizeOverlaps(annot, aln, mode="IntersectionStrict")
##     counts3 <- summarizeOverlaps(annot, aln, mode="IntersectionNotEmpty")

## ----eval=FALSE----------------------------------------------------------
##     synthTrxCountsTable <- data.frame(
##                                  assays(counts1)$counts,
##                                  assays(counts2)$counts,
##                                  assays(counts3)$counts)
##     colnames(synthTrxCountsTable) <- c("union","intStrict","intNotEmpty")
##     rownames(synthTrxCountsTable) <- rownames(counts1)
##     sds <- apply(synthTrxCountsTable,1,sd)
##     sum(sds!=0)
##     sum(sds!=0)/length(sds)
## 
##     synthTrxCountsTable[which.max(sds),]
##     syntheticTrx[which.max(sds),]

## ------------------------------------------------------------------------
    aln <- readGAlignments(
     BamFile(file.path(extdata(),"BAM","207_subset_sortmerna_trimmomatic_sorted.bam")))

## ------------------------------------------------------------------------
    cover <- coverage(aln)

## ----eval=FALSE----------------------------------------------------------
##     cover
##     # this object is compressed to save space. It is an RLE (Running Length Encoding)
##     # we can look at a section of chromosome 4 say between bp 1 and 1000
##     # which gives us the number of read overlapping each of those bases
##     as.vector(cover[["Chr19"]])[6553903:6561936]

## ----eval=FALSE----------------------------------------------------------
##     islands <- slice(cover, 1)
##     islandPeakHeight <- viewMaxs(islands)
##     islandWidth <- width(islands)

## ----eval=FALSE----------------------------------------------------------
##     candidateExons <- islands[islandPeakHeight >= 2L & islandWidth >=150L]
##     candidateExons[["Chr19"]]

## ------------------------------------------------------------------------
    sum(cigarOpTable(cigar(aln))[,"N"] > 0)

## ------------------------------------------------------------------------
    aln[cigarOpTable(cigar(aln))[,"N"] > 0 & seqnames(aln) == "Chr19",]
    aln[cigarOpTable(cigar(aln))[,"N"] > 0 & seqnames(aln) == "Chr19",][3:6,]

## ------------------------------------------------------------------------
    cover[["Chr19"]][131278:132275]

## ------------------------------------------------------------------------
    splice.reads <- aln[cigarOpTable(cigar(aln))[,"N"] > 0 & seqnames(aln) == "Chr19",]
    cherry.pick.sel <- c(13492:13498,13511:13517)
    read.start <- start(splice.reads)[cherry.pick.sel]
    donor.pos <- read.start - 1 +
    as.numeric(sapply(strsplit(cigar(splice.reads)[cherry.pick.sel],"M"),"[",1))
     
    acceptor.pos <- read.start - 1 +
     sapply(
       lapply(
         lapply(strsplit(cigar(splice.reads)[cherry.pick.sel],"M|N"),"[",1:2),
         as.integer),
       sum)

## ------------------------------------------------------------------------
    Chr19 <- readDNAStringSet(file.path(extdata(),
                                       "FASTA",
                                       "Ptrichocarpa_v3.0_210-Chr19.fa.gz"))

## ------------------------------------------------------------------------
    donor <- Views(subject=Chr19[[1]],
                        start=donor.pos-8,
                        end=donor.pos+11)
    acceptor <- Views(subject=Chr19[[1]],
                      start=acceptor.pos-10,
                       end=acceptor.pos+9)

## ------------------------------------------------------------------------
    donor <- DNAStringSet(donor)
    minus.acceptor <- reverseComplement(donor[1:7])

    acceptor <- DNAStringSet(acceptor)
    minus.donor <- reverseComplement(acceptor[1:7])

    donor[1:7] <- minus.donor
    acceptor[1:7] <- minus.acceptor

## ----eval=FALSE----------------------------------------------------------
##     library(seqLogo)
##     pwm <- makePWM(cbind(
##      alphabetByCycle(donor)[c("A","C","G","T"),]/length(donor),
##      alphabetByCycle(acceptor)[c("A","C","G","T"),]/length(acceptor))
##     )
##     seqLogo(pwm)

## ----eval=FALSE----------------------------------------------------------
##     dir(file.path(extdata(),
##                   pattern="^[A,T].*\\.bam$",
##                   full.names=TRUE)

## ----eval=FALSE----------------------------------------------------------
##     system.file(
##        "extdata",
##        "Dmel-mRNA-exon-r5.52.gff3",
##        package="RnaSeqTutorial")

## ----eval=FALSE----------------------------------------------------------
##     ## we start by loading the packages
##     library("easyRNASeq")
##     library("RnaSeqTutorial")
## 
##     ## looking up the function definition
##     show(simpleRNASeq)
## 
##     ## creating the BamFileList
##     bamFileList <- getBamFileList(
##                     dir(path=system.file("extdata",
##                         package="RnaSeqTutorial"),
##                         pattern="^[A,T].*\\.bam$",
##                         full.names=TRUE))
## 
##     ## creating the AnnotParam object
##     annotParam <- AnnotParam(datasource=system.file(
##                            "extdata",
##                            "Dmel-mRNA-exon-r5.52.gff3",
##                            package="RnaSeqTutorial"))
## 
##     ## creating the BamParam object
##     bamParam <- BamParam(paired=FALSE,stranded=FALSE)
## 
##     ## creating the RnaSeqParam
##     rnaSeqParam <- RnaSeqParam(annotParam=annotParam,
##                               bamParam=bamParam,
##                               countBy="exons")
## 
##     ## and we then get a SummarizedExperiment containing the counts table
##     sexp <- simpleRNASeq(
##        bamFiles=bamFileList,
##        param=rnaSeqParam,
##        verbose=TRUE
##        )
## 
##     ## and look at the counts
##     exon.count <- assays(sexp)$exons
## 
##     head(exon.count[rowSums(exon.count)>0,])

## ----eval=FALSE----------------------------------------------------------
##     sexp <- simpleRNASeq(
##        bamFiles=bamFileList,
##        param=RnaSeqParam(annotParam=annotParam),
##        verbose=TRUE
##        )

## ----eval=FALSE----------------------------------------------------------
##     ## the counts
##     head(assays(sexp)[[1]])
##     ## some non empty counts
##     head(assays(sexp)[[1]][rowSums(assays(sexp)[[1]])!=0,])
##     ## the sample info
##     colData(sexp)
##     ## the 'features' info
##     rowRanges(sexp)

## ----eval=FALSE----------------------------------------------------------
##     vignette("easyRNASeq")

