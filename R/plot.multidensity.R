##' Multiple Density plot
##' 
##' Draws multiple density curves within a single plot.
##' 
##' This function calculates the most suitable axis from the 
##' densities calculated from the provided list of numerical 
##' vectors and display all densities in a single plot
##' 
##' @name Multiple density plot
##' @rdname RnaSeqTutorial-plotmultidensity
##' @aliases plot.multidensity plot.multidensity,list-method
##' @param x a list of numeric vectors
##' @param xlab the label of the x axis, default to \code{x}
##' @param col the color palette to use default to \code{brewer.pal(8,"Dark2")}
##' from the \code{\link[RColorBrewer:RColorBrewer]{RColorBrewer}} package.
##' @param legend.x the position of the legend, default to "top", 
##' see the \code{\link[graphics:legend]{legend}} function for 
##' more details about the positional argument
##' @param xlim the xlim defaults to NULL. If NULL the limit is automatically 
##' calculated
##' @param ylim the ylim defaults to NULL. If NULL the limit is automatically 
##' calculated
##' @param lty the line type, default to 1 (blank), see 
##' \code{\link[graphics:par]{par}} for more details 
##' @param legend.cex the character expansion of the legend, defaults to 1, see
##' \code{\link[graphics:par]{cex}} for more details
##' @param legend the legend to be displayed, if NULL and x is a named list, 
##' the names are used as legend, otherwise the legend is omitted. If provided,
##' it must be a list of length \code{length(x)}
##' @param ... additional argument passed to the \code{plot} and \code{lines} 
##' functions.
##' @return NULL
##' @seealso \code{\link[RColorBrewer:RColorBrewer]{RColorBrewer}},
##' \code{\link[graphics:legend]{legend}} and
##' \code{\link[graphics:par]{par}}
##' @examples
##' \dontrun{
##' ## TODO add an object after putting results in the data folder
##' ## and use that for plotting the example
##' data(RST-densities)
##' plot.multidensity(density.list)
##' }
##' 

##' @exportMethod plot.multidensity
setGeneric(name="plot.multidensity",
           def=function(x,xlab="x",col=brewer.pal(8,"Dark2"),
                        legend.x="top",xlim=NULL,ylim=NULL,
                        lty=1,legend.cex=1,legend=NULL,...){
             standardGeneric("plot.multidensity")
           })

setMethod(f="plot.multidensity",
          signature="list",
          definition=function(x,xlab="x",col=brewer.pal(8,"Dark2"),
                              legend.x="top",xlim=NULL,ylim=NULL,
                              lty=1,legend.cex=1,legend=NULL,...){
            
            ## check
            stopifnot(is.list(x))
            
            ## densities
            dens <- lapply(x,density)
            if(is.null(xlim)){
              xlim <- range(sapply(dens,"[[","x"))
            }
            if(is.null(ylim)){
              ylim <- range(sapply(dens,"[[","y"))
            }
            
            ## lty
            if(length(lty)==1){
              lty <- rep(lty,length(dens))
            }
            
            ## plot
            plot(0,0,xlim=xlim,ylim=ylim,ylab="density",type="n",xlab=xlab,...)
            
            ## lines
            sapply(1:length(dens),function(i,dens,col,...){
              lines(dens[[i]],col=col[i],
                    lty=lty[i],...)},dens,col,...)
            
            ## legend
            if(!is.null(legend) | !is.null(names(x))){
              if(is.null(legend)){
                leg <- names(x)    
              } else {
                leg <- legend
              }
              legend(legend.x,col=col[1:length(x)],bty="n",
                     legend=leg,lty=lty,cex=legend.cex)
            }
            return(NULL)
          })
