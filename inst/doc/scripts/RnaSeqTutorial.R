## ----message=FALSE,warning=FALSE,results='hide',echo=FALSE---------------
options(digits=2)

## ----child = "Chapters/01-Convention.Rmd"--------------------------------

## ------------------------------------------------------------------------
a <- 1

## ------------------------------------------------------------------------
(a <- 1)


## ----child = "Chapters/02-Introduction.Rmd"------------------------------



## ----child = "Chapters/03-R-Introduction.Rmd"----------------------------

## ------------------------------------------------------------------------
    ## assign values 5, 4, 3, 2, 1 to variable 'x'
    x <- c(5, 4, 3, 2, 1)
    x

## ------------------------------------------------------------------------
    x[2:4]

## ------------------------------------------------------------------------
    log(x)

## ------------------------------------------------------------------------
    c(1.1, 1.2, 1.3)         # numeric
    c(FALSE, TRUE, FALSE)    # logical
    c("foo", "bar", "baz")   # character, single or double quote ok
    as.character(x)          # convert 'x' to character
    typeof(x)                # the number 5 is numeric, not integer
    typeof(2L)               # append 'L' to force integer
    typeof(2:4)              # ':' produces a sequence of integers

## ------------------------------------------------------------------------
    sex <- factor(c("Male", "Female", NA), levels=c("Female", "Male"))
    sex

## ------------------------------------------------------------------------
    lst <- list(a=1:3, b=c("foo", "bar"), c=sex)
    lst

## ------------------------------------------------------------------------
    lst[c(3, 1)]             # another list
    lst[["a"]]               # the element itself, selected by name

## ------------------------------------------------------------------------
    df <- data.frame(age=c(27L, 32L, 19L),
                     sex=factor(c("Male", "Female", "Male")))
    df
    df[c(1, 3),]
    df[df$age > 20,]

## ------------------------------------------------------------------------
    m <- matrix(1:12, nrow=3)
    m
    m[c(1, 3), c(2, 4)]
    m[, 3]
    m[, 3, drop=FALSE]

## ------------------------------------------------------------------------
    x <- rnorm(1000, sd=1)
    y <- x + rnorm(1000, sd=.5)
    fit <- lm(y ~ x)       # formula describes linear regression 
    fit                    # an 'S3' object
    anova(fit)
    sqrt(var(resid(fit)))  # residuals accessor and subsequent transforms
    class(fit)

## ------------------------------------------------------------------------
    y <- 5:1
    log(y)
    args(log)        # arguments 'x' and 'base'; see ?log
    log(y, base=2)   # 'base' is optional, with default value
    try(log())       # 'x' required; 'try' continues even on error
    args(data.frame) # ... represents variable number of arguments

## ------------------------------------------------------------------------
    log(base=2, y)   # match argument 'base' by name, 'x' by position

## ------------------------------------------------------------------------
    args(anova)
    args(stats:::anova.glm)

## ------------------------------------------------------------------------
      pdataFile <- system.file(package="RnaSeqTutorial", "extdata", "pData.csv")

## ------------------------------------------------------------------------
    pdata <- read.table(pdataFile)  
    dim(pdata)
    names(pdata)
    summary(pdata)

## ------------------------------------------------------------------------
    head(pdata[,"sex"], 3)
    head(pdata$sex, 3)
    head(pdata[["sex"]], 3)
    sapply(pdata, class)

## ------------------------------------------------------------------------
    table(pdata$sex, useNA="ifany")

## ------------------------------------------------------------------------
    with(pdata, table(mol.biol, useNA="ifany"))

## ------------------------------------------------------------------------
    ridx <- pdata$mol.biol %in% c("BCR/ABL", "NEG")

## ------------------------------------------------------------------------
    table(ridx)
    sum(ridx)

## ------------------------------------------------------------------------
    pdata1 <- pdata[ridx,]

## ------------------------------------------------------------------------
    levels(pdata$mol.biol)

## ------------------------------------------------------------------------
    pdata1$mol.biol <- factor(pdata1$mol.biol)
    table(pdata1$mol.biol)

## ------------------------------------------------------------------------
    with(pdata1, t.test(age ~ mol.biol))

## ------------------------------------------------------------------------
    boxplot(age ~ mol.biol, pdata1)

## ------------------------------------------------------------------------
     boxplot(age ~ mol.biol, pdata1,notch=TRUE,ylab="age (yr)",
            main="Age distribution by genotype",xlab="genotype")

## ------------------------------------------------------------------------
    library(lattice)
    plt <- dotplot(variety ~ yield | site, data = barley, groups = year,
                   xlab = "Barley Yield (bushels/acre)" , ylab=NULL,
                   key = simpleKey(levels(barley$year), space = "top", 
                     columns=2),
                   aspect=0.5, layout = c(2,3))
    print(plt)

## ------------------------------------------------------------------------
    length(search())
    search()
    base::log(1:3)

## ----eval=FALSE----------------------------------------------------------
##     library(RnaSeqTutorial)
##     sessionInfo()

## ----eval=FALSE----------------------------------------------------------
##     help.start()

## ----eval=FALSE----------------------------------------------------------
##     ?data.frame
##     ?lm
##     ?anova             # a generic function
##     ?anova.lm          # an S3 method, specialized for 'lm' objects

## ------------------------------------------------------------------------
    methods(anova)
    methods(class="glm")

## ----eval=FALSE----------------------------------------------------------
##     ls
##     getAnywhere("anova.loess")

## ------------------------------------------------------------------------
    utils::head
    methods(head)
    head(head.matrix)

## ----eval=FALSE----------------------------------------------------------
##     library(Biostrings)
##     showMethods(complement)

## ----eval=FALSE----------------------------------------------------------
##     showMethods(class="DNAStringSet", where=search())

## ----eval=FALSE----------------------------------------------------------
##     class ? DNAStringSet
##     method ? "complement,DNAStringSet"

## ----eval=FALSE----------------------------------------------------------
##     selectMethod(complement, "DNAStringSet")

## ----eval=FALSE----------------------------------------------------------
##     vignette(package="RnaSeqTutorial")

## ----eval=FALSE----------------------------------------------------------
##     vignette("RnaSeqTutorial")

## ----eval=FALSE----------------------------------------------------------
##     ?library
##     library(Biostrings)
##     ?alphabetFrequency
##     class?GAlignments
##     vignette(package="GenomicRanges")

## ----eval=FALSE, keep.source=TRUE----------------------------------------
##     # not evaluated
##     colClasses <-
##      c("NULL", "integer", "numeric", "NULL")
##     df <- read.table("myfile", colClasses=colClasses)

## ----keep.source=TRUE----------------------------------------------------
    x <- runif(100000); x2 <- x^2
    m <- matrix(x2, nrow=1000); y <- rowSums(m)

## ----eval=FALSE, keep.source=TRUE----------------------------------------
##     ## not evaluated
##     result <- numeric(nrow(df))
##     for (i in seq_len(nrow(df)))
##      result[[i]] <- some_calc(df[i,])

## ------------------------------------------------------------------------
    unlist(list(a=1:2)) # name 'a' becomes 'a1', 'a2'
    unlist(list(a=1:2), use.names=FALSE)   # no names

## ----eval=FALSE, keep.source=TRUE----------------------------------------
##     ## not evaluated
##     library(limma) # microarray linear models
##     fit <- lmFit(eSet, design)

## ----keep.source=TRUE----------------------------------------------------
    x <- 1:100; s <- sample(x, 10)
    inS <- x %in% s

## ------------------------------------------------------------------------
    m <- matrix(runif(200000), 20000)
    replicate(5, system.time(apply(m, 1, sum))[[1]])
    replicate(5, system.time(rowSums(m))[[1]])

## ------------------------------------------------------------------------
    res1 <- apply(m, 1, sum)
    res2 <- rowSums(m)
    identical(res1, res2)
    identical(c(1, -1), c(x=1, y=-1))
    all.equal(c(1, -1), c(x=1, y=-1),
              check.attributes=FALSE)


## ----child = "Chapters/04-Preprocessing.Rmd"-----------------------------



## ----child = "Chapters/05-Alignment.Rmd"---------------------------------



## ----child = "Chapters/06-Prelude.Rmd"-----------------------------------

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


## ----child = "Chapters/07-R-Alignments.Rmd"------------------------------

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


## ----child = "Chapters/08-R-Annotations.Rmd"-----------------------------

## ----message=FALSE,warning=FALSE,results='hide',echo=FALSE---------------
    options(digits=2)
    library(RnaSeqTutorial)

## ------------------------------------------------------------------------
    library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
    txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene

## ------------------------------------------------------------------------
    set.seed(123)
    fbids <- sample(keys(txdb),10,FALSE,)
    txnm <- select(txdb, fbids, "TXNAME", "GENEID")
    nrow(txnm)
    head(txnm, 3)

## ------------------------------------------------------------------------
    cds <- cdsBy(txdb, "tx", use.names=TRUE)[txnm$TXNAME[6]]
    cds[1]

## ------------------------------------------------------------------------
    library(BSgenome.Dmelanogaster.UCSC.dm3)
    txx <- extractTranscriptSeqs(Dmelanogaster, cds)
    length(txx)
    head(txx, 3)
    head(translate(txx), 3)

## ------------------------------------------------------------------------
library(IRanges)
library(genomeIntervals)

## ------------------------------------------------------------------------
gff <- readGff3(file.path(extdata(),
                          "GFF3/Ptrichocarpa_v3.0_210_gene_exons.gff3.gz"),
                quiet=TRUE)

## ------------------------------------------------------------------------
gff
nrow(gff[gff$type=="exon",])
nrow(gff[gff$type=="mRNA",])
nrow(gff[gff$type=="gene",])

## ------------------------------------------------------------------------
sel <- gff$type == "mRNA"
transcriptGeneMapping <- data.frame(getGffAttribute(gff[sel], "ID"), 
                                    getGffAttribute(gff[sel], "Parent")
)
head(transcriptGeneMapping)

## ------------------------------------------------------------------------
sel <- gff$type=="exon"
rngList<- split(IRanges(start=gff[sel,1],end=gff[sel,2]),
                transcriptGeneMapping[match(
                  sapply(strsplit(getGffAttribute(gff[sel,],"Parent"),","),"[",1),
                  transcriptGeneMapping$ID),"Parent"])
rngList
mostExons <- rev(names(table(elementNROWS(rngList))))[1]
mostExons

## ------------------------------------------------------------------------
rngList<- IRanges::reduce(rngList)
rngList
rev(names(table(elementNROWS(rngList))))[1]

## ------------------------------------------------------------------------
exons <- gff[sel,]
syntheticGeneModel<- exons[rep(
  match(names(rngList),
        transcriptGeneMapping[
          match(sapply(strsplit(getGffAttribute(exons,"Parent"),","),"[",1),
                transcriptGeneMapping$ID),"Parent"]),
  elementNROWS(rngList)),]

## ------------------------------------------------------------------------
syntheticGeneModel[,1]<- unlist(start(rngList))
syntheticGeneModel[,2]<- unlist(end(rngList))
levels(syntheticGeneModel$source)<- "inhouse"

## ------------------------------------------------------------------------
exonNumber<- lapply(elementNROWS(rngList),":",1)
sel<- strand(syntheticGeneModel)[cumsum(elementNROWS(rngList))] == "+"
exonNumber[sel]<- sapply(exonNumber[sel],rev)

## ------------------------------------------------------------------------
syntheticGeneModel$gffAttributes<- paste("ID=",
                                         rep(names(rngList),elementNROWS(rngList)),
                                         ":",unlist(exonNumber),";Parent=",
                                         rep(names(rngList),elementNROWS(rngList)),".0",sep="")

## ------------------------------------------------------------------------
writeGff3(syntheticGeneModel,file="~/Ptrichocarpa_v3.0_210_synthetic_transcripts.gff3")

sel <- syntheticGeneModel$type=="exon"
annot <- split(GRanges(seqnames=seq_name(syntheticGeneModel[sel]),
                       ranges=IRanges(start=syntheticGeneModel[sel,1],
                                      end=syntheticGeneModel[sel,2]),
                       strand=strand(syntheticGeneModel[sel])),
               getGffAttribute(syntheticGeneModel,"Parent")[sel,1]
)

save(annot,file="~/Ptrichocarpa_v3.0_210_synthetic_transcripts.rda")

## ------------------------------------------------------------------------
synthTrx <- createSyntheticTranscripts(
    file.path(extdata(),"GFF3/Ptrichocarpa_v3.0_210_gene_exons.gff3.gz"),
    verbose=FALSE)


## ----child = "Chapters/09-R-Interlude.Rmd"-------------------------------

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


## ----child = "Chapters/10-Counting.Rmd"----------------------------------



## ----child = "Chapters/11-R-Counts.Rmd"----------------------------------

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


## ----child = "Chapters/12-R-Biological-QA.Rmd"---------------------------

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


## ----child = "Chapters/13-R-Differential-Expression.Rmd"-----------------

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


## ----child = "Commons/17-SessionInfo.Rmd"--------------------------------

## ----session info, echo=FALSE--------------------------------------------
sessionInfo()


## ----child = "Commons/14-Acknowledgments.Rmd"----------------------------



## ----child = "Commons/15-Footnotes.Rmd"----------------------------------



## ----child = "Commons/16-Images.Rmd"-------------------------------------



