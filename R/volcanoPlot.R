##' Volcano plot
##' 
##' Draws a Volcano plot, i.e. the differential expression log 
##' odds vs. expression log fold change between condition. 
##' 
##' The function has been devised to draw a MA plot from a \code{DESeq}
##' result object.
##' 
##' @name Volcano plot
##' @rdname RnaSeqTutorial-volcanoPlot
##' @aliases volcanoPlot volcanoPlot,DataFrame-method
##' @param object an object of class \code{\linkS4class{DataFrame}} having at
##' least 2 columns named \code{log2FoldChange} and \code{padj}
##' @param alpha the significance cutoff to be used when plotting, defaults to
##' 10e-3.
##' @return NULL
##' @seealso \code{\link[DESeq:newCountDataSet]{DESeq package}}
##' @examples
##' \dontrun{
##' ## TODO add an object after putting results in the data folder
##' ## and use that for plotting the example
##' data(RST-DataFrame)
##' volcanoPlot(dataFrame)
##' }
##' 

##' @exportMethod volcanoPlot
setGeneric(name="volcanoPlot",def=function(object,alpha=0.001){
  standardGeneric("volcanoPlot")
})

setMethod(f="volcanoPlot",
          signature="DataFrame",
          definition=function(object,alpha=0.001){
            
            ## selectors
            sel <- ! is.na(object$padj)
            sel2 <- object$padj[sel]<=alpha 
            
            ## plot
            heatscatter(object$log2FoldChange[sel],
                        -log10(object$padj[sel]),
                        main="Volcano",xlab="Log2 Fold Change", 
                        ylab="- log(10) adj. p-value")
            
            ## legend
            legend("topleft",bty="n",paste("cutoff @",alpha),lty=2,col="gray")
            
            ## points
            points(object$log2FoldChange[sel][sel2],
                   -log10(object$padj[sel][sel2]),col="lightblue",pch=19)
            points(object$log2FoldChange[sel][sel2],
                   -log10(object$padj[sel][sel2]),col="dodgerblue3",
                   pch=19,cex=0.5)
            
            ## circle the points for the dot plot
            abline(h=-log10(alpha),lty=2,col="gray")
            
            return(NULL)
          })
