## ------------------------------------------------------------------------
    ## assign values 5, 4, 3, 2, 1 to variable 'x'
    x <- c(5, 4, 3, 2, 1)
    x

## ------------------------------------------------------------------------
    x[2:4]

## ------------------------------------------------------------------------
    log(x)

## ------------------------------------------------------------------------
    c(1.1, 1.2, 1.3)         # numeric
    c(FALSE, TRUE, FALSE)    # logical
    c("foo", "bar", "baz")   # character, single or double quote ok
    as.character(x)          # convert 'x' to character
    typeof(x)                # the number 5 is numeric, not integer
    typeof(2L)               # append 'L' to force integer
    typeof(2:4)              # ':' produces a sequence of integers

## ------------------------------------------------------------------------
    sex <- factor(c("Male", "Female", NA), levels=c("Female", "Male"))
    sex

## ------------------------------------------------------------------------
    lst <- list(a=1:3, b=c("foo", "bar"), c=sex)
    lst

## ------------------------------------------------------------------------
    lst[c(3, 1)]             # another list
    lst[["a"]]               # the element itself, selected by name

## ------------------------------------------------------------------------
    df <- data.frame(age=c(27L, 32L, 19L),
                     sex=factor(c("Male", "Female", "Male")))
    df
    df[c(1, 3),]
    df[df$age > 20,]

## ------------------------------------------------------------------------
    m <- matrix(1:12, nrow=3)
    m
    m[c(1, 3), c(2, 4)]
    m[, 3]
    m[, 3, drop=FALSE]

## ------------------------------------------------------------------------
    x <- rnorm(1000, sd=1)
    y <- x + rnorm(1000, sd=.5)
    fit <- lm(y ~ x)       # formula describes linear regression 
    fit                    # an 'S3' object
    anova(fit)
    sqrt(var(resid(fit)))  # residuals accessor and subsequent transforms
    class(fit)

## ------------------------------------------------------------------------
    y <- 5:1
    log(y)
    args(log)        # arguments 'x' and 'base'; see ?log
    log(y, base=2)   # 'base' is optional, with default value
    try(log())       # 'x' required; 'try' continues even on error
    args(data.frame) # ... represents variable number of arguments

## ------------------------------------------------------------------------
    log(base=2, y)   # match argument 'base' by name, 'x' by position

## ------------------------------------------------------------------------
    args(anova)
    args(stats:::anova.glm)

## ------------------------------------------------------------------------
      pdataFile <- system.file(package="RnaSeqTutorial", "extdata", "pData.csv")

## ------------------------------------------------------------------------
    pdata <- read.table(pdataFile)  
    dim(pdata)
    names(pdata)
    summary(pdata)

## ------------------------------------------------------------------------
    head(pdata[,"sex"], 3)
    head(pdata$sex, 3)
    head(pdata[["sex"]], 3)
    sapply(pdata, class)

## ------------------------------------------------------------------------
    table(pdata$sex, useNA="ifany")

## ------------------------------------------------------------------------
    with(pdata, table(mol.biol, useNA="ifany"))

## ------------------------------------------------------------------------
    ridx <- pdata$mol.biol %in% c("BCR/ABL", "NEG")

## ------------------------------------------------------------------------
    table(ridx)
    sum(ridx)

## ------------------------------------------------------------------------
    pdata1 <- pdata[ridx,]

## ------------------------------------------------------------------------
    levels(pdata$mol.biol)

## ------------------------------------------------------------------------
    pdata1$mol.biol <- factor(pdata1$mol.biol)
    table(pdata1$mol.biol)

## ------------------------------------------------------------------------
    with(pdata1, t.test(age ~ mol.biol))

## ------------------------------------------------------------------------
    boxplot(age ~ mol.biol, pdata1)

## ------------------------------------------------------------------------
     boxplot(age ~ mol.biol, pdata1,notch=TRUE,ylab="age (yr)",
            main="Age distribution by genotype",xlab="genotype")

## ------------------------------------------------------------------------
    library(lattice)
    plt <- dotplot(variety ~ yield | site, data = barley, groups = year,
                   xlab = "Barley Yield (bushels/acre)" , ylab=NULL,
                   key = simpleKey(levels(barley$year), space = "top", 
                     columns=2),
                   aspect=0.5, layout = c(2,3))
    print(plt)

## ------------------------------------------------------------------------
    length(search())
    search()
    base::log(1:3)

## ----eval=FALSE----------------------------------------------------------
##     library(RnaSeqTutorial)
##     sessionInfo()

## ----eval=FALSE----------------------------------------------------------
##     help.start()

## ----eval=FALSE----------------------------------------------------------
##     ?data.frame
##     ?lm
##     ?anova             # a generic function
##     ?anova.lm          # an S3 method, specialized for 'lm' objects

## ------------------------------------------------------------------------
    methods(anova)
    methods(class="glm")

## ----eval=FALSE----------------------------------------------------------
##     ls
##     getAnywhere("anova.loess")

## ------------------------------------------------------------------------
    utils::head
    methods(head)
    head(head.matrix)

## ----eval=FALSE----------------------------------------------------------
##     library(Biostrings)
##     showMethods(complement)

## ----eval=FALSE----------------------------------------------------------
##     showMethods(class="DNAStringSet", where=search())

## ----eval=FALSE----------------------------------------------------------
##     class ? DNAStringSet
##     method ? "complement,DNAStringSet"

## ----eval=FALSE----------------------------------------------------------
##     selectMethod(complement, "DNAStringSet")

## ----eval=FALSE----------------------------------------------------------
##     vignette(package="RnaSeqTutorial")

## ----eval=FALSE----------------------------------------------------------
##     vignette("RnaSeqTutorial")

## ----eval=FALSE----------------------------------------------------------
##     ?library
##     library(Biostrings)
##     ?alphabetFrequency
##     class?GAlignments
##     vignette(package="GenomicRanges")

## ----eval=FALSE, keep.source=TRUE----------------------------------------
##     # not evaluated
##     colClasses <-
##      c("NULL", "integer", "numeric", "NULL")
##     df <- read.table("myfile", colClasses=colClasses)

## ----keep.source=TRUE----------------------------------------------------
    x <- runif(100000); x2 <- x^2
    m <- matrix(x2, nrow=1000); y <- rowSums(m)

## ----eval=FALSE, keep.source=TRUE----------------------------------------
##     ## not evaluated
##     result <- numeric(nrow(df))
##     for (i in seq_len(nrow(df)))
##      result[[i]] <- some_calc(df[i,])

## ------------------------------------------------------------------------
    unlist(list(a=1:2)) # name 'a' becomes 'a1', 'a2'
    unlist(list(a=1:2), use.names=FALSE)   # no names

## ----eval=FALSE, keep.source=TRUE----------------------------------------
##     ## not evaluated
##     library(limma) # microarray linear models
##     fit <- lmFit(eSet, design)

## ----keep.source=TRUE----------------------------------------------------
    x <- 1:100; s <- sample(x, 10)
    inS <- x %in% s

## ------------------------------------------------------------------------
    m <- matrix(runif(200000), 20000)
    replicate(5, system.time(apply(m, 1, sum))[[1]])
    replicate(5, system.time(rowSums(m))[[1]])

## ------------------------------------------------------------------------
    res1 <- apply(m, 1, sum)
    res2 <- rowSums(m)
    identical(res1, res2)
    identical(c(1, -1), c(x=1, y=-1))
    all.equal(c(1, -1), c(x=1, y=-1),
              check.attributes=FALSE)

