##' Data accessor for the RnaSeqTutorial package
##' 
##' Simply lookup the extdata path.
##' 
##' This function simply lookup the extdata path in the installation directory
##' of the RnaSeqTutorial package
##' 
##' @name Package data accessor
##' @rdname RnaSeqTutorial-extdata
##' @aliases extdata extdata,missing-method
##' @param x missing
##' @return The file path to the RnaSeqTutorial package extdata folder
##' @examples
##'   extdata()
##'   

#' @exportMethod extdata
setGeneric(name="extdata",
           def=function(x){
             standardGeneric("extdata")
           })

setMethod(f="extdata",
          signature="missing",
          definition=function(x){
            system.file("extdata",package="RnaSeqTutorial")
          })

