## to re-build the documentation and namespace
## library(roxygen2)
## roxygenize("RnaSeqTutorial")

###==========================
### package details
###==========================
##' Count summarization and normalization pipeline for Next Generation
##' Sequencing data.
##'
##' Offers functionalities to summarize read counts per feature of interest,
##' e.g. exons, transcripts, genes, etc.
##' Offers functionalities to normalize the summarized counts using 3rd party
##' packages like \code{\link[DESeq:newCountDataSet]{DESeq}}
##' or \code{\link[edgeR:DGEList]{edgeR}}.
##'
##'
##' @name RnaSeqTutorial package
##' @rdname RnaSeqTutorial-package
##' @aliases RnaSeqTutorial-package
##' @docType package
##' @author Nicolas Delhomme, Bastian Schiffthaler, Ismael Padioleau
##' @keywords package
##' @seealso
##'     \code{\link[easyRNASeq:easyRNASeq-package]{easyRNASeq}}
##'
##' @examples
##'     \dontrun{
##'     library("RnaSeqTutorial")
##'   TODO put some examples there
##'     }
##'
NULL

###==========================
## To define the NAMESPACE
###==========================
## library(codetoolsBioC)
## library(RnaSeqTutorial)
## writeNamespaceImports("RnaSeqTutorial")
## import while namespaces
##' @import BSgenome.Dmelanogaster.UCSC.dm3 GenomicAlignments
##' GenomicFeatures genomeIntervals easyRNASeq Rsamtools
##' TxDb.Dmelanogaster.UCSC.dm3.ensGene leeBamViews
## import classes
##' @importClassesFrom DESeq CountDataSet
##' @importClassesFrom methods numeric
##' @importClassesFrom S4Vectors DataFrame
## import S4 methods
## @importMethodsFrom Biobase fData varMetadata
##' @importMethodsFrom BiocGenerics counts density paste plotMA subset
## @importMethodsFrom BSgenome colnames nrow
## @importMethodsFrom easyRNASeq plotDispLSD
##' @importMethodsFrom IRanges pmax pmin quantile
##' @importMethodsFrom S4Vectors colnames nrow
##' @importMethodsFrom seqLogo plot
##' @importMethodsFrom ShortRead "%in%" lapply sapply
## import methods
##' @importFrom DESeq fitInfo
##' @importFrom graphics abline axis axTicks box contour grid image legend lines
##' mtext par points
##' @importFrom grDevices topo.colors
##' @importFrom LSD heatscatter
##' @importFrom RColorBrewer brewer.pal
##' @importFrom utils tail
## and export
## @exportClass BamFileList
## @exportMethod assay
## @export chromosomeFilter
NULL

## To check
## library(BiocCheck)
## BiocCheck("../RnaSeqTutorial")
