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

## ----fqc, eval=FALSE, engine="bash"--------------------------------------
##   # First, make an output directory called "qa" and a sub-directory of
##   # that called "raw"
##   mkdir -p qa/raw
## 
##   # We can also make a directory to hold the raw data
##   mkdir raw
## 
##   # Then make copies of your sample files into your home folder. Always try to
##   # work on copies of the raw source data.
##   cp share/Day1/data/fastq/*.fq.gz raw
## 
##   # Then you can run fastqc on the file(s) you have selected. FastQC can run
##   # a number of files in parallel (see the -t option)
##   fastqc -o qa/raw -t 16 raw/*.fq.gz
## 
##   # Finally, also run multiqc to summarize all data in one plot
##   cd qa/raw
##   multiqc .

## ---- eval=FALSE, engine="bash"------------------------------------------
##     # First uncompress your forward and reverse reads:
##     gunzip read_1.fq.gz read_2.fq.gz
## 
##     # Merge forward and reverse reads into a single file (this is a requirement
##     # of the sortmerna software)
##     merge-paired-reads.sh read_1.fq read_2.fq read-interleaved.fq
## 
##     # Run sortmerna
##     sortmerna --ref $SORTMERNA_DB --reads read-interleaved.fq --aligned \
##     read-rRNA-hits --other read-sortmerna --log -a 16 -v --paired_in --fastx
## 
##     # Split the file again in two separate files, one for forward, one for reverse
##     # orientation
##     unmerge-paired-reads.sh read-sortmerna.fq read-sortmerna_1.fq read-sortmerna_2.fq
## 
##     # Delete the interleaved files and compress the result files to save space
##     rm read-interleaved.fq read-sortmerna.fq
##     gzip read-sortmerna_1.fq read-sortmerna_2.fq

## ---- eval=FALSE, engine="bash"------------------------------------------
##     #!/bin/bash
## 
##     READ_FW="$1"
##     READ_RV="$2"
## 
##     FILEBASE=$(basename "${READ_FW/_1.fq.gz/}")
## 
##     echo "Uncompressing FASTQ data of $FILEBASE"
##     gunzip "$READ_FW" "$READ_RV"
## 
##     READ_FW="${READ_FW%.gz}"
##     READ_RV="${READ_RV%.gz}"
## 
##     echo "Merging pairs of $FILEBASE"
##     merge-paired-reads.sh "$READ_FW" "$READ_RV" "${FILEBASE}_interleaved.fq"
## 
##     echo "Running SortMeRNA for $FILEBASE"
##     sortmerna --ref $SORTMERNA_DB --reads "${FILEBASE}_interleaved.fq" --aligned \
##     "${FILEBASE}-rRNA-hits" --other  "${FILEBASE}-sortmerna" --log -a 16 \
##     -v --paired_in --fastx
## 
##     echo "Unmerging SortMeRNA filtered pairs for $FILEBASE"
##     unmerge-paired-reads.sh "${FILEBASE}-sortmerna.fq" \
##     "${FILEBASE}-sortmerna_1.fq" "${FILEBASE}-sortmerna_2.fq"
## 
##     echo "Doing cleanup for $FILEBASE"
##     gzip "$READ_FW" "$READ_RV" "${FILEBASE}-sortmerna_1.fq" \
##     "${FILEBASE}-sortmerna_2.fq" "${FILEBASE}-rRNA-hits.fq"
##     rm "${FILEBASE}_interleaved.fq" "${FILEBASE}-sortmerna.fq"

## ---- eval=FALSE, engine="bash"------------------------------------------
##     mkdir ~/sortmerna
##     cd ~/sortmerna
##     find ../raw -name "*.fq.gz"  | sort | head -n 4 | while read READ_FW
##     do
##       read READ_RV
##       bash ../runSortMeRNA.sh $READ_FW $READ_RV
##     done

## ---- eval=FALSE, engine="bash"------------------------------------------
##    cp ~/share/Day1/data/sortmerna/* .

## ---- eval=FALSE, engine="bash"------------------------------------------
##     cd ~/sortmerna
##     multiqc .

## ---- eval=FALSE, engine="bash"------------------------------------------
##     mkdir ~/qa/sortmerna
##     fastqc -o ~/qa/sortmerna -t 16 ~/sortmerna/*sortmerna*.fq.gz

## ---- eval=FALSE, engine="bash"------------------------------------------
##     mkdir ~/trimmomatic
##     cd ~/trimmomatic
##     find ../sortmerna -name "*sortmerna_[12].fq.gz" | sort | head -n 4 | while read FW_READ
##     do
##       read RV_READ
##       FILEBASE=$(basename "${FW_READ%_1.fq.gz}")
##       trimmomatic PE -threads 16 -phred64 "$FW_READ" "$RV_READ" \
##       "$FILEBASE-trimmomatic_1.fq.gz" "$FILEBASE-trimmomatic-unpaired_1.fq.gz" \
##       "$FILEBASE-trimmomatic_2.fq.gz" "$FILEBASE-trimmomatic-unpaired_2.fq.gz"  \
##       ILLUMINACLIP:"/usr/share/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa":2:30:10 \
##       SLIDINGWINDOW:5:20 MINLEN:50 2> "$FILEBASE-timmomatic.log"
##     done

## ---- eval=FALSE, engine="bash"------------------------------------------
##     cd ~/trimmomatic
##     cp ~/share/Day1/data/trimmomatic/* .
##     multiqc .

## ---- eval=FALSE, engine="bash"------------------------------------------
##     mkdir ~/qa/trimmomatic
##     fastqc -o ~/qa/trimmomatic -t 16 ~/trimmomatic/*trimmomatic_[12].fq.gz

## ---- eval=FALSE, engine="bash"------------------------------------------
##     mkdir ~/qc_report
##     ln -s ~/qa/raw/*zip ~/qc_report
##     ln -s ~/qa/sortmerna/*zip ~/qc_report
##     ln -s ~/qa/trimmomatic/*zip ~/qc_report
##     ln -s ~/sortmerna/*.log ~/qc_report
##     ln -s ~/trimmomatic/*.log ~/qc_report


## ----child = "Chapters/05-Alignment.Rmd"---------------------------------

## ---- eval=FALSE, engine="bash"------------------------------------------
##     mkdir ~/star
##     cd ~/star
##     find ../trimmomatic -name "*trimmomatic_[12].fq.gz" | sort | head -n 4 | while read FW_READ
##     do
##       read RV_READ
##       FILEBASE=$(basename "${FW_READ%_1.fq.gz}")
##       STAR --genomeDir ~/share/Day1/data/indices/STAR --readFilesIn "$FW_READ" "$RV_READ" \
##         --runThreadN 16 --alignIntronMax 11000 --outSAMstrandField intronMotif \
##         --sjdbGTFfile ~/share/Day1/data/reference/gff/Ptrichocarpa_v3.0_210_synthetic-gene-models-wo-introns.gff3 \
##         --readFilesCommand zcat --outFileNamePrefix "$FILEBASE-STAR"  --outSAMmapqUnique 254 \
##         --outQSconversionAdd -31 --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate \
##         --quantMode TranscriptomeSAM --outFilterMultimapNmax 100 --chimSegmentMin 1 --outWigType bedGraph \
##         --sjdbGTFtagExonParentTranscript Parent
##     done

## ---- eval=FALSE, engine="bash"------------------------------------------
##     samtools view ~/star/202_subset-sortmerna-trimmomatic-STARAligned.sortedByCoord.out.bam | head

## ---- eval=FALSE, engine="bash"------------------------------------------
##     samtools index ~/star/202_subset-sortmerna-trimmomatic-STARAligned.sortedByCoord.out.bam

## ---- eval=FALSE, engine="bash"------------------------------------------
##     mkdir ~/star_logs
##     cp ~/star/*Log.final.out ~/star_logs
##     cd ~/star_logs
##     multiqc .

## ---- eval=FALSE, engine="bash"------------------------------------------
##     cd
##     kallisto index -i Potri03-synthetic-transcripts.idx \
##     ~/share/Day1/data/reference/fasta/Potri03-synthetic-transcripts.fa

## ---- eval=FALSE, engine="bash"------------------------------------------
##     mkdir ~/kallisto
##     cd ~/kallisto
##     find ../trimmomatic -name "*trimmomatic_[12].fq.gz" | sort | head -n 4 | while read FW_READ
##     do
##       read RV_READ
##       FILEBASE=$(basename "${FW_READ%_1.fq.gz}")
##       kallisto quant -i ../Potri03-synthetic-transcripts.idx -b 100 \
##       -o . -t 16 "$FW_READ" "$RV_READ"
##       # Kallisto doesn't let us specify an output filename so we rename all
##       # output files
##       mv "abundance.tsv" $FILEBASE-"abundance.tsv"
##       mv "abundance.h5" $FILEBASE-"abundance.h5"
##       mv "run_info.json" $FILEBASE-"run_info.json"
##     done


## ----child = "Chapters/07-Linux-Alignments.Rmd"--------------------------

## ---- eval=FALSE, engine="bash"------------------------------------------
## mkdir ~/alignments
## cp ~/share/Day1/data/star/{202,207,213.1,221,226.1}_subset*Aligned.sortedByCoord.out.bam ~/alignments

## ---- eval=FALSE, engine="bash"------------------------------------------
## cd ~/alignments
## samtools view 202_subset-sortmerna-trimmomatic-STARAligned.sortedByCoord.out.bam | less -S

## ---- eval=FALSE, engine="bash"------------------------------------------
## samtools

## ---- eval=FALSE, engine="bash"------------------------------------------
## samtools sort <BAM> > <BAM_sorted>

## ---- eval=FALSE, engine="bash"------------------------------------------
## find . -name "*.bam" | xargs -I{} bash -c 'samtools sort $0 > ${0/.bam/_sorted.bam}' {}

## ---- eval=FALSE, engine="bash"------------------------------------------
## samtools index <BAM>

## ---- eval=FALSE, engine="bash"------------------------------------------
## find . -name "*sorted.bam" | xargs -I {} samtools index {}

## ---- eval=FALSE, engine="bash"------------------------------------------
## samtools flagstat <BAM>

## ---- eval=FALSE, engine="bash"------------------------------------------
## samtools stats <BAM> | less -S

## ---- eval=FALSE, engine="bash"------------------------------------------
## samtools idxstats <BAM> | less -S

## ---- eval=FALSE, engine="bash"------------------------------------------
## find . -name "*sorted.bam" | xargs -I{} bash -c 'samtools flagstat $0 > ${0/.bam/_flagstat.txt}' {}
## 
## find . -name "*sorted.bam" | xargs -I{} bash -c 'samtools idxstats $0 > ${0/.bam/_idxstats.txt}' {}
## 
## find . -name "*sorted.bam" | xargs -I{} bash -c 'samtools stats $0 > ${0/.bam/_stats.txt}' {}
## 
## multiqc .

## ---- eval=FALSE, engine="bash"------------------------------------------
## bedtools

## ---- eval = FALSE, engine="bash"----------------------------------------
## mkdir ~/bedtools
## cd ~/bedtools
## bedtools multicov -bams ../star/*.sortedByCoord.out.bam -bed ~/share/Day1/data/bed/putative_sex_locus.bed | less -S


## ----child = "Chapters/08-Linux-R-Annotations.Rmd"-----------------------

## ----message=FALSE,warning=FALSE,results='hide',echo=FALSE---------------
    options(digits=2)
    library(RnaSeqTutorial)

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
    nrow(unique(gff[gff$type=="exon",])) / nrow(gff[gff$type=="gene",])

## ------------------------------------------------------------------------
sel <- duplicated(gff[gff$type=="exon",])
firstPosition <- which(sel)[1]
firstDuplicate <- gff[gff$type=="exon",][firstPosition,]

## ------------------------------------------------------------------------
sel <- gff$type=="exon"
gff[sel,][seqnames(gff[sel,]) == seqnames(firstDuplicate) &
              gff[sel,1] == firstDuplicate[,1] &
              gff[sel,2] == firstDuplicate[,2],]

## ------------------------------------------------------------------------
synthTrx <- createSyntheticTranscripts(
    file.path(extdata(),"GFF3/Ptrichocarpa_v3.0_210_gene_exons.gff3.gz"),
    verbose=FALSE)

## ---- eval=FALSE---------------------------------------------------------
## writeGff3(synthTrx,file="~/gff3/Ptrichocarpa_v3.0_210_synthetic_transcripts.gff3")

## ---- eval=FALSE---------------------------------------------------------
## library(DESeq2)
## setwd("~/share/CIBNOR2017/Day1/data/htseq-comp/")
## tab <- data.frame("SampleID"=sub("_subset.txt","",dir("full-gff/")),
##                     "File"=dir("full-gff/"),
##                     "Sample"=sub("_subset.txt","",dir("full-gff/")))
## 
## full <- DESeqDataSetFromHTSeqCount(tab,"full-gff",design = ~Sample)
## 
## synth <- DESeqDataSetFromHTSeqCount(tab,"synth-gff",design = ~Sample)
## 
## par(mfrow=c(1,2))
## barplot(colSums(counts(synth)),ylim=c(1,1e6),las=2,main="Synthetic Transcripts gff3")
## barplot(colSums(counts(full)),ylim=c(1,1e6),las=2,main="Original Transcripts gff3")

## ------------------------------------------------------------------------
    library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
    txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene

## ------------------------------------------------------------------------
    set.seed(123)
    fbids <- sample(keys(txdb),10,FALSE)
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
synthTrx<- exons[rep(
  match(names(rngList),
        transcriptGeneMapping[
          match(sapply(strsplit(getGffAttribute(exons,"Parent"),","),"[",1),
                transcriptGeneMapping$ID),"Parent"]),
  elementNROWS(rngList)),]

## ------------------------------------------------------------------------
synthTrx[,1]<- unlist(start(rngList))
synthTrx[,2]<- unlist(end(rngList))
levels(synthTrx$source)<- "inhouse"

## ------------------------------------------------------------------------
exonNumber<- lapply(elementNROWS(rngList),":",1)
sel<- strand(synthTrx)[cumsum(elementNROWS(rngList))] == "+"
exonNumber[sel]<- sapply(exonNumber[sel],rev)

## ------------------------------------------------------------------------
synthTrx$gffAttributes<- paste("ID=",
                                         rep(names(rngList),elementNROWS(rngList)),
                                         ":",unlist(exonNumber),";Parent=",
                                         rep(names(rngList),elementNROWS(rngList)),".0",sep="")


## ----child = "Chapters/10-Counting.Rmd"----------------------------------



## ----child = "Chapters/12-R-Biological-QA.Rmd"---------------------------

## ----message=FALSE,warning=FALSE,results='hide'--------------------------
    library(RnaSeqTutorial)

## ------------------------------------------------------------------------
res <- mclapply(dir(file.path(extdata(),"htseq"),
                    pattern="^[2,3].*_STAR\\.txt",
                    full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=16)
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
sex <- samples$sex[match(names(res),samples$sample)]
date <- factor(samples$date[match(names(res),samples$sample)])

## ------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(
                       sex=sex,
                       date=date),
  design = ~ date + sex)

## ------------------------------------------------------------------------
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- names(res)
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

## ---- eval=FALSE---------------------------------------------------------
##     library(tximport)
##     files <- list.files("~/kallisto",
##                     pattern = ".*-abundance.tsv",
##                     full.names = TRUE)
##     tx <- suppressMessages(tximport(files = files,
##                                      type = "kallisto",
##                                      txOut = TRUE))

## ---- eval=FALSE---------------------------------------------------------
## k.countTable <- round(tx$counts)

## ---- eval=FALSE---------------------------------------------------------
## library(LSD)
## par(mfrow=c(3,3))
## dev.null <- sapply(1:ncol(count.table),function(i){
##     heatscatter(log2(count.table[,i]+1),
##                 log2(k.countTable[,i]+1),
##                 main=colnames(count.table)[i],
##                 xlab="HTSeq count (log2 + pseudo-count)",
##                 ylab="Kallisto count (log2 + pseudo-count)",
##                 )
## })


## ----child = "Commons/17-SessionInfo.Rmd"--------------------------------

## ----session info, echo=FALSE--------------------------------------------
sessionInfo()


## ----child = "Commons/14-Acknowledgments.Rmd"----------------------------



## ----child = "Commons/15-Footnotes.Rmd"----------------------------------



## ----child = "Commons/16-Images.Rmd"-------------------------------------



