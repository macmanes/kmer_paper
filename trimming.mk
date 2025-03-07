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
TRIM=0
CPU=2
BCPU=$(CPU)
RUN=test123
READ1=left.fq
READ2=right.fq
BCODES=barcodes.fa
TRIMMOMATIC= $(shell which trimmomatic-0.32.jar)
JELLY= $(shell which jellyfish)
FASTOOL= $(shell which fastool)
TRINITY= $(shell which Trinity)
TRANS=$(shell locate FL_trans_analysis_pipeline.pl)
MUS=$(shell find /share -name Mus_musculus.GRCm38.71.cdna.all.fa 2> /dev/null)
PFAM=$(shell find /share -name  Pfam-AB.hmm.bin 2> /dev/null)
CONVERT=$(shell locate fq2fa.py)
LIMIT=60
Z=5

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}

.PHONY: check clean
all: check $(RUN)_left.fastq $(RUN)_right.fastq $(RUN)_right.fastq.goodkmers.fa \
	$(RUN)_left.fastq.goodkmers.fa $(RUN)_right.fastq.kmerfilt $(RUN)_left.fastq.kmerfilt $(READ2).kmerfilt.fa \
	$(READ1).kmerfilt.fa $(RUN)/jellyfish.kmers.fa $(RUN)/both.fa $(RUN)/Trinity.fasta \
	$(RUN).Trinity.fasta.pslx $(RUN).Trinity.fasta.pep




$(RUN)_left.fastq $(RUN)_right.fastq: $(READ1) $(READ2)
	@echo TIMESTAMP: `date +'%a %d%b%Y  %H:%M:%S'` About to start trimming
		mkdir $(RUN)
		java -Xmx$(MEM)g -jar $(TRIMMOMATIC) PE -baseout $(RUN).fq -phred$(PHRED) -threads $(CPU) \
		$(READ1) \
		$(READ2) \
		ILLUMINACLIP:${MAKEDIR}/$(BCODES):2:40:15 \
		LEADING:$(TRIM) TRAILING:$(TRIM) SLIDINGWINDOW:4:$(TRIM) MINLEN:$(MINLEN) | tee trim.log ;
		cat $(RUN)_1P.fq $(RUN)_1U.fq > $(RUN)_left.fastq ;
		cat $(RUN)_2P.fq $(RUN)_2U.fq > $(RUN)_right.fastq ;
		rm $(RUN)_1P.fq $(RUN)_1U.fq $(RUN)_2P.fq $(RUN)_2U.fq ;
		@echo TIMESTAMP: `date +'%a %d%b%Y  %H:%M:%S'` Finished trimming '\n\n'


$(RUN)_right.fastq.goodkmers.fa $(RUN)_left.fastq.goodkmers.fa $(RUN)_right.fastq.kmerfilt $(RUN)_left.fastq.kmerfilt:$(RUN)_left.fastq $(RUN)_right.fastq
	@echo TIMESTAMP: `date +'%a %d%b%Y  %H:%M:%S'` location aware trimming
	python /share/khmer/scripts/extract-untrusted-kmers.py --quiet -k 25 --limit $(LIMIT) -Z $(Z) $(RUN)_left.fastq &
	python /share/khmer/scripts/extract-untrusted-kmers.py --quiet -k 25 --limit $(LIMIT) -Z $(Z) $(RUN)_right.fastq
	@echo TIMESTAMP: `date +'%a %d%b%Y  %H:%M:%S'` ending location aware trimming

$(READ2).kmerfilt.fa $(READ1).kmerfilt.fa:$(RUN)_right.fastq.kmerfilt $(RUN)_left.fastq.kmerfilt
	python $(CONVERT) $(RUN)_left.fastq.kmerfilt $(READ1).kmerfilt.fa &
	python $(CONVERT) $(RUN)_right.fastq.kmerfilt $(READ2).kmerfilt.fa 

$(RUN)/jellyfish.kmers.fa:$(RUN)_right.fastq.goodkmers.fa $(RUN)_left.fastq.goodkmers.fa $(READ1).kmerfilt.fa $(READ2).kmerfilt.fa
	@echo TIMESTAMP: `date +'%a %d%b%Y  %H:%M:%S'` starting jellyfish
	$(JELLY) count -m25 -t$(CPU) -C -s3G -o jf.out $(RUN)_right.fastq.goodkmers.fa $(RUN)_left.fastq.goodkmers.fa $(READ1).kmerfilt.fa $(READ2).kmerfilt.fa
	#$(JELLY) merge jf.out* -o merged.jf 
	$(JELLY) dump -o $(RUN)/jellyfish.kmers.fa jf.out
	rm jf.out $(RUN)_right.fastq.goodkmers.fa $(RUN)_left.fastq.goodkmers.fa $(READ1).kmerfilt.fa $(READ2).kmerfilt.fa &
	@echo TIMESTAMP: `date +'%a %d%b%Y  %H:%M:%S'` ending jellyfish

$(RUN)/both.fa:$(RUN)_left.fastq $(RUN)_right.fastq
	@echo TIMESTAMP: `date +'%a %d%b%Y  %H:%M:%S'` starting fastool
	$(FASTOOL) --illumina-trinity --to-fasta $(RUN)_left.fastq >> left.fa
	$(FASTOOL) --illumina-trinity --to-fasta $(RUN)_right.fastq >> right.fa
	cat left.fa right.fa > $(RUN)/both.fa
	rm left.fa right.fa &
	touch $(RUN)/jellyfish.1.finished

$(RUN)/Trinity.fasta:$(RUN)_left.fastq $(RUN)_right.fastq
	@echo TIMESTAMP: `date +'%a %d%b%Y  %H:%M:%S'` starting trinity
	$(TRINITY) --min_kmer_cov $(MINK) --seqType $(SEQ) --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G --bflyCPU $(BCPU) \
	--left $(RUN)_left.fastq --right $(RUN)_right.fastq --group_pairs_distance 999 --CPU $(CPU) --output $(RUN) | tee $(RUN).trinity.pe.log

$(RUN).Trinity.fasta.pslx:$(RUN)/Trinity.fasta
	$(TRANS) --target $(MUS) --query $(RUN)/Trinity.fasta; rm *maps *selected *summary
$(RUN).Trinity.fasta.pep:$(RUN)/Trinity.fasta
	TransDecoder --CPU $(CPU) -t $(RUN)/Trinity.fasta \
	--search_pfam $(PFAM) | tee pfam10.log; \
	rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep $(RUN).Trinity.fasta.pep

check:
	@echo TIMESTAMP: `date +'%a %d%b%Y  %H:%M:%S'` ---Begin--- '\n\n'
	@echo "\n\n\n"###I am checking to see if you have all the dependancies installed.### "\n"
	command -v bwa mem >/dev/null 2>&1 || { echo >&2 "I require BWA but it's not installed.  Aborting."; exit 1; }
	@echo BWA is Installed
	command -v samtools >/dev/null 2>&1 || { echo >&2 "I require samTools but it's not installed.  Aborting."; exit 1; }
	@echo samTools is Installed
	command -v $(TRINITY) >/dev/null 2>&1 || { echo >&2 "I require Trinity but it's not installed.  Aborting."; exit 1; }
	@echo Trinity is Installed
	if [ -f $(READ1) ]; then echo 'left fastQ exists'; else echo 'Im having trouble finding your left fastQ file, check PATH \n'; exit 1; fi;
	if [ -f $(READ2) ]; then echo 'right fastQ exists \n'; else echo 'Im having trouble finding your right fastQ file, check PATH \n'; fi;
	chmod -w $(READ1) 2>/dev/null; true
	chmod -w $(READ2) 2>/dev/null; true
