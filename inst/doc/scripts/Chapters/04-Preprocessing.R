## ----fqc, eval=FALSE, engine="bash"--------------------------------------
##   # First, make an output directory called "qa" and a sub-directory of
##   # that called "raw"
##   mkdir -p qa/raw
## 
##   # We can also make a directory to hold the raw data
##   mkdir raw
## 
##   # Then make copies of your sample files into your home folder. Always try to
##   # work on copies of the raw source data.
##   cp share/Day1/data/fastq/*.fq.gz raw
## 
##   # Then you can run fastqc on the file(s) you have selected. FastQC can run
##   # a number of files in parallel (see the -t option)
##   fastqc -o qa/raw -t 16 raw/*.fq.gz
## 
##   # Finally, also run multiqc to summarize all data in one plot
##   cd qa/raw
##   multiqc .

## ---- eval=FALSE, engine="bash"------------------------------------------
##     # First uncompress your forward and reverse reads:
##     gunzip read_1.fq.gz read_2.fq.gz
## 
##     # Merge forward and reverse reads into a single file (this is a requirement
##     # of the sortmerna software)
##     merge-paired-reads.sh read_1.fq read_2.fq read-interleaved.fq
## 
##     # Run sortmerna
##     sortmerna --ref $SORTMERNA_DB --reads read-interleaved.fq --aligned \
##     read-rRNA-hits --other read-sortmerna --log -a 16 -v --paired_in --fastx
## 
##     # Split the file again in two separate files, one for forward, one for reverse
##     # orientation
##     unmerge-paired-reads.sh read-sortmerna.fq read-sortmerna_1.fq read-sortmerna_2.fq
## 
##     # Delete the interleaved files and compress the result files to save space
##     rm read-interleaved.fq read-sortmerna.fq
##     gzip read-sortmerna_1.fq read-sortmerna_2.fq

## ---- eval=FALSE, engine="bash"------------------------------------------
##     #!/bin/bash
## 
##     READ_FW="$1"
##     READ_RV="$2"
## 
##     FILEBASE=$(basename "${READ_FW/_1.fq.gz/}")
## 
##     echo "Uncompressing FASTQ data of $FILEBASE"
##     gunzip "$READ_FW" "$READ_RV"
## 
##     READ_FW="${READ_FW%.gz}"
##     READ_RV="${READ_RV%.gz}"
## 
##     echo "Merging pairs of $FILEBASE"
##     merge-paired-reads.sh "$READ_FW" "$READ_RV" "${FILEBASE}_interleaved.fq"
## 
##     echo "Running SortMeRNA for $FILEBASE"
##     sortmerna --ref $SORTMERNA_DB --reads "${FILEBASE}_interleaved.fq" --aligned \
##     "${FILEBASE}-rRNA-hits" --other  "${FILEBASE}-sortmerna" --log -a 16 \
##     -v --paired_in --fastx
## 
##     echo "Unmerging SortMeRNA filtered pairs for $FILEBASE"
##     unmerge-paired-reads.sh "${FILEBASE}-sortmerna.fq" \
##     "${FILEBASE}-sortmerna_1.fq" "${FILEBASE}-sortmerna_2.fq"
## 
##     echo "Doing cleanup for $FILEBASE"
##     gzip "$READ_FW" "$READ_RV" "${FILEBASE}-sortmerna_1.fq" \
##     "${FILEBASE}-sortmerna_2.fq" "${FILEBASE}-rRNA-hits.fq"
##     rm "${FILEBASE}_interleaved.fq" "${FILEBASE}-sortmerna.fq"

## ---- eval=FALSE, engine="bash"------------------------------------------
##     mkdir ~/sortmerna
##     cd ~/sortmerna
##     find ../raw -name "*.fq.gz"  | sort | head -n 4 | while read READ_FW
##     do
##       read READ_RV
##       bash ../runSortMeRNA.sh $READ_FW $READ_RV
##     done

## ---- eval=FALSE, engine="bash"------------------------------------------
##    cp ~/share/Day1/data/sortmerna/* .

## ---- eval=FALSE, engine="bash"------------------------------------------
##     cd ~/sortmerna
##     multiqc .

## ---- eval=FALSE, engine="bash"------------------------------------------
##     mkdir ~/qa/sortmerna
##     fastqc -o ~/qa/sortmerna -t 16 ~/sortmerna/*sortmerna*.fq.gz

## ---- eval=FALSE, engine="bash"------------------------------------------
##     mkdir ~/trimmomatic
##     cd ~/trimmomatic
##     find ../sortmerna -name "*sortmerna_[12].fq.gz" | sort | head -n 4 | while read FW_READ
##     do
##       read RV_READ
##       FILEBASE=$(basename "${FW_READ%_1.fq.gz}")
##       trimmomatic PE -threads 16 -phred64 "$FW_READ" "$RV_READ" \
##       "$FILEBASE-trimmomatic_1.fq.gz" "$FILEBASE-trimmomatic-unpaired_1.fq.gz" \
##       "$FILEBASE-trimmomatic_2.fq.gz" "$FILEBASE-trimmomatic-unpaired_2.fq.gz"  \
##       ILLUMINACLIP:"/usr/share/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa":2:30:10 \
##       SLIDINGWINDOW:5:20 MINLEN:50 2> "$FILEBASE-timmomatic.log"
##     done

## ---- eval=FALSE, engine="bash"------------------------------------------
##     cd ~/trimmomatic
##     cp ~/share/Day1/data/trimmomatic/* .
##     multiqc .

## ---- eval=FALSE, engine="bash"------------------------------------------
##     mkdir ~/qa/trimmomatic
##     fastqc -o ~/qa/trimmomatic -t 16 ~/trimmomatic/*trimmomatic_[12].fq.gz

## ---- eval=FALSE, engine="bash"------------------------------------------
##     mkdir ~/qc_report
##     ln -s ~/qa/raw/*zip ~/qc_report
##     ln -s ~/qa/sortmerna/*zip ~/qc_report
##     ln -s ~/qa/trimmomatic/*zip ~/qc_report
##     ln -s ~/sortmerna/*.log ~/qc_report
##     ln -s ~/trimmomatic/*.log ~/qc_report

