##' Density plot adapted from the LSD package
##' 
##' Draws a topographic density plot.
##' 
##' This function was extracted from the LSD package and
##' was modified to suppress the change of margins
##' 
##' @name Density plot
##' @rdname RnaSeqTutorial-densityPlot
##' @aliases densityPlot densityPlot,numeric,numeric-method
##' @param x a numeric vector for the x axis value
##' @param y a numeric vector for the y axis value
##' @param grid the size of the grid, default to 100
##' @param ncol the number of color to use, default to 30 
##' @param nlevels the number of levels passed as argument to the 
##' \code{\link[graphics:contour]{contour}} function, defaults to 10 
##' @param ... additional argument passed to the 
##' \code{\link[graphics:image]{image}} function.
##' @return NULL
##' @seealso \code{\link[graphics:contour]{contour}},
##' \code{\link[graphics:image]{image}} and
##' \code{\link[LSD:heatscatter]{heatscatter}}
##' @examples
##' \dontrun{
##' ## TODO add an object after putting results in the data folder
##' ## and use that for plotting the example
##' data(RST-DataFrame)
##' sel <- ! is.na(dataFrame$padj)
##' densityPlot(log10(dataFrame$baseMean),
##'             dataFrame$log2FoldChange)
##' }

##' @exportMethod densityPlot 
setGeneric(name="densityPlot",
           def=function(x, y, grid = 100, 
                        ncol = 30, 
                        nlevels = 10, ...){
             standardGeneric("densityPlot")
           })

setMethod(f="densityPlot",
          signature=c("numeric","numeric"),
          definition=function(x, y, grid = 100, ncol = 30, 
                              nlevels = 10, ...){
            if (!is.vector(x) | !is.vector(y)) 
              stop("First two argument must be vectors !")
            if (length(x) != length(y)) 
              stop("Data vectors must be of the same length !")
            
            ## function
            "kde2d.adj" <- function(x, y, h, n = 25, 
                                    lims = c(range(x), 
                                             range(y)), only = "none") {
              nx = length(x)
              gx = seq.int(lims[1], lims[2], length.out = n)
              gy = seq.int(lims[3], lims[4], length.out = n)
              bandwidth.nrd.adj = function(x) {
                r = quantile(x, c(0.25, 0.75))
                h = (r[2] - r[1])/1.34
                return(4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5))
              }
              if (missing(h)) {
                bx = bandwidth.nrd.adj(x)
                by = bandwidth.nrd.adj(y)
                if (all(c(bx, by) == 0)) {
                  h = rep(0.01, 2)
                }
                else if (any(c(bx, by) == 0)) {
                  h = rep(max(bx, by), 2)
                }
                else {
                  h = c(bx, by)
                }
              }
              else h = rep(h, length.out = 2)
              h = h/4
              ax = outer(gx, x, "-")/h[1]
              ay = outer(gy, y, "-")/h[2]
              norm.ax = dnorm(ax)
              norm.ay = dnorm(ay)
              if (only == "x") {
                norm.ay = rep(1, length(ay))
              }
              if (only == "y") {
                norm.ax = rep(1, length(ax))
              }
              z = tcrossprod(matrix(norm.ax, , nx), 
                             matrix(norm.ay, , nx))/(nx * h[1] * h[2])
              list(x = gx, y = gy, z = z)
            }
            
            ## main
            d = kde2d.adj(x, y, n = grid)
            z <- d$z
            nrz <- nrow(z)
            ncz <- ncol(z)
            couleurs <- tail(topo.colors(trunc(1.4 * ncol)), ncol)
            image(d, col = couleurs, ...)
            contour(d, add = TRUE, nlevels = nlevels)
            box()
            return(NULL)
          })
