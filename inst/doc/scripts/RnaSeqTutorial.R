## ----message=FALSE,warning=FALSE,results='hide',echo=FALSE---------------
options(digits=2)

## ----child = "Chapters/01-Convention.Rmd"--------------------------------

## ------------------------------------------------------------------------
a <- 1

## ------------------------------------------------------------------------
(a <- 1)


## ----child = "Chapters/02-Introduction.Rmd"------------------------------



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
## cp ~/share/Day1/star/{202,207,213.1,221,226.1}_subset* ~/alignments

## ---- eval=FALSE, engine="bash"------------------------------------------
## cd ~/alignments
## samtools view 202_subset_sortmerna_trimmomatic_sorted.bam | less -S

## ---- eval=FALSE, engine="bash"------------------------------------------
## samtools

## ---- eval=FALSE, engine="bash"------------------------------------------
## samtools sort <BAM> > <BAM_sorted>

## ---- eval=FALSE, engine="bash"------------------------------------------
## find . -name "*.bam" | xargs -I{} bash -c 'samtools sort $0 > ${0/.bam/_sorted.bam}' {}

## ---- eval=FALSE, engine="bash"------------------------------------------
## samtools index <BAM>

## ---- eval=FALSE, engine="bash"------------------------------------------
## find . -name "*sorted.bam" | xargs -I{} samtools index {}

## ---- eval=FALSE, engine="bash"------------------------------------------
## samtools flagstat <BAM>

## ---- eval=FALSE, engine="bash"------------------------------------------
## samtools stats <BAM> | less -S

## ---- eval=FALSE, engine="bash"------------------------------------------
## samtools idxstats <BAM> | less -S

## ---- eval=FALSE, engine="bash"------------------------------------------
## find . -name "*sorted.bam" | xargs -I{} samtools index {}
## 
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

## ---- eval=FALSE, engine="bash"------------------------------------------
## cd ~/bedtools
## bedtools subtract -a ../star/202_subset-sortmerna-trimmomatic-STARAligned.sortedByCoord.out.bam -b ~/share/Day1/data/reference/gff/Ptrichocarpa_v3.0_210_synthetic-gene-models-wo-introns.gff3 > not_in_feature.bam
## samtools view not_in_feature.bam | less -S


## ----child = "Chapters/08-Linux-R-Annotations.Rmd"-----------------------

## ----b6801, eval=FALSE, engine="bash"------------------------------------
##     cut -f3 Day1/... | sort | uniq -c

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

## ---- eval=FALSE, engine="bash"------------------------------------------
##     cd && mkdir gff3
##     gt gff3 -force -tidy yes -addids yes -fixregionboundaries yes \
##     -retainids yes -sort yes -checkids yes \
##     -o Ptrichocarpa_v3.0_210_gene_exons-validated.gff3 \
##     share/Day1/data/reference/gff/Ptrichocarpa_v3.0_210_gene_exons.gff3.gz \
##     2>&1  | grep -v "##sequence-region"

## ------------------------------------------------------------------------
synthTrx <- createSyntheticTranscripts(
    file.path(extdata(),"GFF3/Ptrichocarpa_v3.0_210_gene_exons.gff3.gz"),
    verbose=FALSE)

## ---- eval=FALSE---------------------------------------------------------
## writeGff3(synthTrx,file="~/gff3/Ptrichocarpa_v3.0_210_synthetic_transcripts.gff3")

## ---- eval=FALSE, engine="bash"------------------------------------------
##     cd && mkdir gff3
##     gt gff3 -force -tidy yes -addids yes -fixregionboundaries yes \
##     -retainids yes -sort yes -checkids yes \
##     -o Ptrichocarpa_v3.0_210_synthetic_transcripts-validated.gff3 \
##     share/Day1/data/reference/gff/Ptrichocarpa_v3.0_210_synthetic_transcripts.gff3 \
##     2>&1  | grep -v "##sequence-region"

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



## ----child = "Commons/17-SessionInfo.Rmd"--------------------------------

## ----session info, echo=FALSE--------------------------------------------
sessionInfo()


## ----child = "Commons/14-Acknowledgments.Rmd"----------------------------



## ----child = "Commons/15-Footnotes.Rmd"----------------------------------



## ----child = "Commons/16-Images.Rmd"-------------------------------------



