# Makefile to use knitr for package vignettes
# None of the rendering engine allow us to use the questMasteR

## The Sources
## -----------
# The placeholder
PDFS = RnaSeqTutorial.pdf 

# The chapters
# for portability we would need to list the files instead
CHAPTERS = $(wildcard Chapters/*.Rmd)

## The Targets
## -----------
.PHONY: usage

# The placeholder; easy
%.tex: %.Rnw
	@echo "Sweave $<"
	@cd `dirname $<`; "$(R_HOME)/bin/R" CMD Sweave `basename $<`; cd ..

%.pdf: %.tex
	@echo "PDF $<"
	$(R_HOME)/bin/Rscript -e "tools::texi2pdf('$*.tex')"

## The chapters
# We want them incorporated in the template, the template knitted
# the output copied to inst/doc/Chapters. Next, we want the chapters 
# trimmed of the inv blocks and stored in tmp/

%.html: %.Rmd
	@echo "Knitting chapter $<"
	@cp $< ../inst/doc/Chapters
	@sed "s:TOKEN:$*.Rmd:" Templates/ChapterTemplate.Rmd > ../inst/doc/Tmp.Rmd
	@cd ../inst/doc; $(R_HOME)/bin/Rscript -e "rmarkdown::render('Tmp.Rmd')"; mv Tmp.html $(<:%.Rmd=%.html)
	@rm ../inst/doc/Tmp.Rmd
	@cp $< tmp/tmp.Rmd
	@cd tmp; perl ../scripts/printoutCleanup.pl < tmp.Rmd > `basename $<`; rm tmp.Rmd

# The R script
# We purl the tmp Rmd files and output the R script stripped of the invisible
# blocs in inst/doc/scripts/Chapters
%.R: %.Rmd
	@echo "Purling chapter $<"
	@cd tmp; $(R_HOME)/bin/Rscript -e "knitr::purl('`basename $<`')"
	@mv tmp/`basename $(<:%.Rmd=%.R)` ../inst/doc/scripts/Chapters

## Global
all: setup $(PDFS) $(CHAPTERS:%.Rmd=%.html) $(CHAPTERS:%.Rmd=%.R) tutorial clean

setup:
	@echo "Setting Up"
	@mkdir -p tmp
	@mkdir -p ../inst/doc/Chapters
	@rm -f ../inst/doc/Chapters/*.html
	@mkdir -p ../inst/doc/scripts/Chapters
	@rm -f ../inst/doc/scripts/Chapters/*.R
	@cp -r Bibliography/References.bib Commons Images css js header.html questMaster.R ../inst/doc

clean:
	@rm -rf *.tex *.bbl *.blg *.aux *.out *.log *.spl *.md
	@cd ../inst/doc; rm -rf References.bib Commons Images css js header.html questMaster.R
	@cd ../inst/doc/Chapters; rm -f *.Rmd

pristine: clean
	@rm -rf RnaSeqTutorial.pdf

usage:
	@if [ -z $(R_HOME) ]; then \
	  echo "usage:\n \
		$(R_HOME)/bin/R CMD make \n \
		SWEAVE_OPTIONS=\"eval=FALSE\" $(R_HOME)/bin/R CMD make\n \
		make clean\n \
		make pristine" ; \
	  exit 1 ; \
	fi

# The vignette in html
# Use the trimmed chapters to build
tutorial: 
	@mv Chapters Chapters.orig
	@mv tmp Chapters
	@cp Templates/Tutorial.Rmd .
	$(R_HOME)/bin/Rscript -e "rmarkdown::render('Tutorial.Rmd')"
	$(R_HOME)/bin/Rscript -e "knitr::purl('Tutorial.Rmd')"
	@rm -rf Chapters Tutorial.Rmd
	@mv Chapters.orig Chapters
	@mv Tutorial.html ../inst/doc/RnaSeqTutorial.html
	@mv Tutorial.R ../inst/doc/scripts/RnaSeqTutorial.R
