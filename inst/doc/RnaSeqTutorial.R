### R code from vignette source 'RnaSeqTutorial.Rnw'

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: setup
###################################################
stopifnot(BiocInstaller:::BIOC_VERSION == "3.1")
## check for unix as some dependencies are unix dependent
## or would be hard to get working on windows machines
stopifnot(.Platform$OS.type=="unix")


###################################################
### code chunk number 3: open the URL (eval = FALSE)
###################################################
## browseURL(system.file(package="RnaSeqTutorial","doc","RnaSeqTutorial.html"))


