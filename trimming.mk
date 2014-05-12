#!/usr/bin/make -rRsf

###########################################
###        -usage 'assembly.mk RUN=run CPU=2 MEM=15 READ1=/location/of/read1.fastq READ2=/location/of/read2.fastq'
###         -RUN= name of run
###
###
###         -Make sure your samTools, BWA and Trinity are installed and in
###          your path
###
###
############################################


##### No more Editing should be necessary below this line  #####

MINLEN=25
PHRED=33
SEQ=fq
MINK=1
MEM=2
TRIM=2
CPU=2
BCPU=$(CPU)
RUN=test123
READ1=left.fq
READ2=right.fq
BCODES=barcodes.fa
TRIMMOMATIC= $(shell which trimmomatic-0.32.jar)

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}

.PHONY: check clean
all: check $(RUN)_left.$(TRIM).fastq $(RUN)_right.$(TRIM).fastq $(RUN).Trinity.fasta $(RUN).xprs
trim: $(READ1).goodkmers.fa $(READ2).fastq.goodkmers.fa
jelly: check $(RUN).Trinity.fasta
express: check $(RUN).xprs


$(READ1).goodkmers.fa $(READ2).fastq.goodkmers.fa:$(READ1) $(READ2)
	python ~/khmer/scripts/extract-untrusted-kmers.py --quiet -k 25 --limit 60 -Z 5  $(READ1) &
	python ~/khmer/scripts/extract-untrusted-kmers.py --quiet -k 25 --limit 60 -Z 5  $(READ2)

