\name{RnaSeqTutorial}
\alias{RnaSeqTutorial}
\alias{vignette}
\docType{data}
\title{An RNA-Seq tutorial}
\description{
This tutorial, described in the attached RnaSeqTutorial vignette, introduces
the core Bioconductor packages necessary for processing RNA sequencing data.
It then introduces the easyRNAseq package that can simplify this task and 
offer different processing to be done one Next-Generation Sequencing data. 
It was used for various Bioconductor workshops since June 2010.
}
\arguments{
	\item{RnaSeqTutorial}{The vignette containing the tutorial.}
}
\format{
	The different data available are: 
	\itemize{
		\item{in data:}
		\itemize{
			\item \code{gAnnot} A file containing the Drosophila melanogaster genic annotation retrieved from FlyBase v5.29 (June 2010) and converted into a \link[IRanges:RangedData]{RangedData} object. Although the corresponding gff file has been deprecated - see below - that annotation is still usable.
			\item \code{gRngList} A file containing the Drosophila melanogaster mRNA and exon annotation retrieved from FlyBase v5.52 (April 2013) and converted into a \link[GenomicRanges:GRangesList]{GRangesList} object.
		}

		\item{in extdata:}
		\itemize{
			\item \code{ACACTG.bam, ACACTG.bam, ATGGCT.bam, TTGCGA.bam} 4 RNA-Seq samples from D. melanogaster demultiplexed, with their associated .bai indexes
			\item \code{Dmel-mRNA-exon-r5.52.gff3} The gff3 file containing the information stored in the gRngList \link[GenomicRanges:GRangesList]{GRangesList} object mentioned above.
			\item \code{annot.gff} Deprecated. The gff file containing the information stored in the gAnnot \link[IRanges:RangedData]{RangedData} object mentioned above.
			\item \code{gapped.bam} A D.melanogaster RNA-Seq file containing a few examples of gapped alignment (done by tophat) surrounding the Mef2 gene locus
			\item \code{multiplex_export.txt.gz} The original export file containing the multiplexed data mentioned above.
			\item \code{subset_export.txt.gz} A D.melanogaster RNA-Seq export file, restricted to 100,000 reads
			\item \code{subset.bam} The same aligned against the reference genome (D.melanogaster v3).
		}
	}
}
\value{
  Some fantastic knowledge.
}
\examples{
	\dontrun{vignette("RnaSeqTutorial")}
}
\keyword{datasets}

