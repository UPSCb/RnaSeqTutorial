##' Heatscatter based plot
##' 
##' A based on the \code{\link[LSD:demotour]{LSD}} package
##' \code{\link[LSD:heatscatter]{heatscatter}} function. It is an
##' adaptation of the MA plot from the \code{DESeq} package 
##' \code{\link[DESeq:plotMA]{plotMA}} function.
##' 
##' The plot is a modification that integrate a third dimension as a colouring
##' shade from blue (sparse) to yellow (dense). Grey represent outlier data 
##' points. This is based on the \code{\link[LSD:heatscatter]{heatscatter}} 
##' function of the \code{\link[LSD:demotour]{LSD}} package.
##'
##' @name Heatscatter based plot
##' @rdname RnaSeqTutorial-heatscatter
##' @aliases plotMALSD plotMALSD,data.frame-method
##' @param cex the character extension, default to 0.45
##' @param col the color to mark the points exceeding the significant threshold.
##' Default to 'forestgreen'
##' @param linecol the line color, default to"#00000080"
##' @param log whether to report the axis on a log scale. default to "x" and 
##' "xy" for the MA and dispersion plot, respectively. See 
##' \code{\link[graphics:par]{par}} for more details.
##' @param sign the significance threshold cutoff (on the adjusted p-values). 
##' Default to 5 percent
##' @param x a data.frame object, the result of the \code{DESeq} analysis. 
##' It must contain the following columns: \code{baseMean}, 
##' \code{log2FoldChange} and \code{padj}.
##' @param xlab the x axis label. Default to "mean of normalized counts".
##' @param ylab the y axis lable. Default to "dispersion" or 
##' "expression(log[2] ~ fold ~ change)" for the dispertion and MA plot, 
##' respectively.
##' @param ylim the y axis limits
##' @param ... additional argument passed to the 
##' \code{\link[LSD:heatscatter]{heatscatter}} function.
##' @return NULL
##' @seealso \code{\link[LSD:demotour]{LSD}},
##' \code{\link[LSD:heatscatter]{heatscatter}},
##' \code{\link[DESeq:plotDispEsts]{plotDispEsts}},
##' \code{\link[DESeq:plotMA]{plotMA}}
##' @examples
##' \dontrun{
##' ## TODO add an object after putting results in the data folder
##' ## and use that for plotting the example
##' 
##' data(RST-DataFrame)
##' plotMALSD(dataFrame)
##' }
##' 

##' @exportMethod plotMALSD
setGeneric(name="plotMALSD",
           def=function(
             x, ylim, sign=0.05, 
             col = 'forestgreen',
             linecol = "#00000080", 
             xlab = "mean of normalized counts",
             ylab = expression(log[2] ~ fold ~ change), 
             log = "x", cex = 0.45,...){
             standardGeneric("plotMALSD")
           })

setMethod(f="plotMALSD",
          signature="data.frame",
          definition=function(
            x, ylim, sign=0.05, col = 'forestgreen',
            linecol = "#00000080", xlab = "mean of normalized counts",
            ylab = expression(log[2] ~ fold ~ change), log = "x", cex = 0.45, ...){
            if (!(is.data.frame(x) && all(c("baseMean", "log2FoldChange") %in%
                                            colnames(x))))
              stop("'x' must be a data frame with columns named 'baseMean',
             'log2FoldChange'.")
            x <- x[x$baseMean != 0,]
            py = x$log2FoldChange
            if (missing(ylim))
              ylim = c(-1, 1) * quantile(abs(py[is.finite(py)]), probs = 0.99) *
              1.1
            heatscatter(log10(x$baseMean), pmax(ylim[1], pmin(ylim[2], py)),
                        pch = ifelse(py < ylim[1], 6, ifelse(py > ylim[2], 2, 16)),
                        cexplot = cex,
                        xlab = xlab,
                        ylab = ylab,
                        xaxt = 'n',
                        ylim = ylim, ...)
            
            # Fix logged x-axis
            atx <- axTicks(1)
            labels <- sapply(atx, function (i)
              as.expression(bquote(10^ .(i)))
            )
            axis(1, at=atx, labels=labels)
            abline(h = 0, lwd = 4, col = linecol, lty=1)
            
            # Mark the significant DEGs, quite ugly, right?
            sign.df <- x[x$padj <= sign,]
            pointy <- sign.df$log2FoldChange
            points(log10(sign.df$baseMean),
                   pmax(ylim[1], pmin(ylim[2], pointy)),
                   pch = 1,
                   col = col,
                   cex = cex + 0.5)
            return(NULL)
          })
