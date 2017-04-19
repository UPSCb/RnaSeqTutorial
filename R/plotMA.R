##' MA (median vs. average) plot
##' 
##' Draws an MA plot, in the present case (the log2 fold changes vs. the log10 
##' of the mean)
##' 
##' This function has been devised to draw a series of MA plot from a 
##' \code{DESeq} result object.
##' 
##' @name MA plot
##' @rdname RnaSeqTutorial-plotMA
##' @aliases plotMA plotMA,DataFrame-method
##' @param object an object of class \code{\linkS4class{DataFrame}} having at
##' least 3 columns named \code{baseMean}, \code{log2FoldChange} and \code{padj}
##' @param alpha the significance cutoff to be used when plotting, defaults to
##' 10e-3.
##' @return TRUE invisibly
##' @seealso \code{\link[DESeq:newCountDataSet]{DESeq package}}
##' @examples
##' \dontrun{
##' ## TODO add an object after putting results in the data folder
##' ## and use that for plotting the example
##' data(RST-DataFrame)
##' volcanoPlot(dataFrame)
##' }
##' 

##' @exportMethod plotMA
setMethod(f="plotMA",
          signature="DataFrame",
          definition=function(object,alpha=0.001){
            
            ## selectors
            sel <- ! is.na(object$padj)
            sel2 <- object$padj[sel]<=alpha 
            
            ## graphic params
            orig.par <- par(no.readonly=TRUE)
            par(mfrow=c(2,1))
            
            ## plots
            densityPlot(log10(object$baseMean[sel]),
                        object$log2FoldChange[sel],
                        grid=250,ncol=30,nlevels=10,
                        main="MA density estimation"
            )
            mtext("log10 mean expression",side=2,line=2)
            mtext("log2 FC",side=1,line=2)
            
            heatscatter(log10(object$baseMean[sel]),
                        object$log2FoldChange[sel],
                        add.contour=TRUE,main="MA",
                        xlab="log10 mean expression",
                        ylab="log2 FC",sub=paste(sum(sel2),
                                                 "sign. feats. @",
                                                 alpha,"cutoff"))
            
            points(log10(object$baseMean[sel][sel2]),
                   object$log2FoldChange[sel][sel2],
                   col="darkred",pch=19,cex=.5)
            
            legend("topright",pch=19,col="darkred","sign. feats.")
            
            par(orig.par,no.readonly=TRUE)
            invisible(TRUE)
          })
