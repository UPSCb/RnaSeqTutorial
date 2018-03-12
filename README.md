# Description
A tutorial R package that describes the pre-processing of RNA-Seq
data - from the raw (FASTQ) data to the alignment (mapping) to a reference
genome - and its subsequent differential expression analysis at the gene
and transcript level.

# Installation

# Preparation
Go to the vignettes folder

1. In the Templates folder:
    1. In Tutorial.Rmd: modify the title and decide which chapters should be included
    2. In ChapterTemplate.Rmd: modify the title
2. In the Chapters foolder:
    1. Modify/Add Rmd files   
    2. Modify the path of the data (_e.g._ DayXX to match the course schedule)
3. In the vignettes folder, run
    ```{bash, eval=FALSE}
        R CMD make
    ```
If it works, save the html page you need, as well as the Rscript from the _RnaSeqTutorial/inst/doc_ directory.

# Caveats
inv block only numbers inv-XXXX
When building if it fails at the stage at which it creates the RnaSeqTutorial.html
page, be careful to restore the original structure of the vignette folder (_i.e._ copy
Chapters.orig back to Chapters)