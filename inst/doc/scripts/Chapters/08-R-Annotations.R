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

