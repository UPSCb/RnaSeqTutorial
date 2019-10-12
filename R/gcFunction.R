##' GC calculation
##' 
##' Simple function to calculate the GC percentage of sequences
##' 
##' This function simply lookup the amount of G and C bp in DNAStringSet and
##' return the sequences GC proportion. The function was originally written
##' by Martin Morgan for the original version of that package.
##' 
##' @name GC calculation
##' @rdname RnaSeqTutorial-gcFunction
##' @aliases gcFunction gcFunction,DNAStringSet-method
##' @param x a DNAStringSet
##' @return The percent GC of the sequences in \code{x}
##' @examples
##' \dontrun{
##'   gcFunction(someDNAStringSet)
##' }

#' @exportMethod gcFunction
setGeneric(name="gcFunction",
           def=function(x){
             standardGeneric("gcFunction")
           })

setMethod(f="gcFunction",
          signature="DNAStringSet",
          definition=function(x){
            alf <- alphabetFrequency(x, as.prob=TRUE)
            rowSums(alf[,c("G", "C")])
          })
