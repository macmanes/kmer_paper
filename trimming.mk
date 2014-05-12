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
LIMIT=60
Z=5

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}

.PHONY: check clean
all: check $(RUN)_left.fastq $(RUN)_right.fastq $(RUN)_right.fastq.goodkmers.fa \
	$(RUN)_left.fastq.goodkmers.fa $(RUN)_right.fastq.kmerfilt $(RUN)_left.fastq.kmerfilt $(READ2).kmerfilt.fa \
	$(READ1).kmerfilt.fa $(RUN)/jellyfish.kmers.fa $(RUN)/both.fa $(RUN).Trinity.fasta



#trim: $(RUN)_left.$(TRIM).fastq $(RUN)_right.$(TRIM).fastq
#khmer:$(READ1).goodkmers.fa $(READ2).fastq.goodkmers.fa
#jelly: check $(RUN).Trinity.fasta
#express: check $(RUN).xprs



$(RUN)_left.fastq $(RUN)_right.fastq: $(READ1) $(READ2)
	@echo TIMESTAMP: `date +'%a %d%b%Y  %H:%M:%S'` About to start trimming
		mkdir $(RUN)
		java -Xmx$(MEM)g -jar $(TRIMMOMATIC) PE -baseout $(RUN).fq -phred$(PHRED) -threads $(CPU) \
		$(READ1) \
		$(READ2) \
		ILLUMINACLIP:${MAKEDIR}/$(BCODES):2:40:15 \
		LEADING:$(TRIM) TRAILING:$(TRIM) SLIDINGWINDOW:4:$(TRIM) MINLEN:$(MINLEN) 2> trim.log ;
		cat $(RUN)_1P.fq $(RUN)_1U.fq > $(RUN)_left.fastq ;
		cat $(RUN)_2P.fq $(RUN)_2U.fq > $(RUN)_right.fastq ;
		rm $(RUN)_1P.fq $(RUN)_1U.fq $(RUN)_2P.fq $(RUN)_2U.fq ;
		@echo TIMESTAMP: `date +'%a %d%b%Y  %H:%M:%S'` Finished trimming '\n\n'


$(RUN)_right.fastq.goodkmers.fa $(RUN)_left.fastq.goodkmers.fa $(RUN)_right.fastq.kmerfilt $(RUN)_left.fastq.kmerfilt:$(RUN)_left.fastq $(RUN)_right.fastq
	python ~/khmer/scripts/extract-untrusted-kmers.py --quiet -k 25 --limit $(LIMIT) -Z $(Z) $(RUN)_left.fastq &
	python ~/khmer/scripts/extract-untrusted-kmers.py --quiet -k 25 --limit $(LIMIT) -Z $(Z) $(RUN)_right.fastq


$(READ2).kmerfilt.fa $(READ1).kmerfilt.fa:$(RUN)_right.fastq.kmerfilt $(RUN)_left.fastq.kmerfilt
	python ~/Desktop/python/fq2fa.py $(READ1).kmerfilt $(READ1).kmerfilt.fa &
	python ~/Desktop/python/fq2fa.py $(READ2).kmerfilt $(READ2).kmerfilt.fa 

$(RUN)/jellyfish.kmers.fa:$(READ1).goodkmers.fa $(READ2).fastq.goodkmers.fa $(READ1).kmerfilt.fa $(READ2).kmerfilt.fa
	$(JELLY) count -m25 -t&(CPU) -C -s3G -o jf.out $(READ1).goodkmers.fa $(READ2).fastq.goodkmers.fa $(READ1).kmerfilt.fa $(READ2).kmerfilt.fa
	$(JELLY) dump -o $(RUN)/jellyfish.kmers.fa merged.jf

$(RUN)/both.fa:$(RUN)_left.fastq $(RUN)_right.fastq
	$(FASTOOL) --rev  --illumina-trinity --to-fasta $(RUN)_left.fastq >> left.fa
	$(FASTOOL) --illumina-trinity --to-fasta $(RUN)_right.fastq >> right.fa
	cat left.fa right.fa > $(RUN)/both.fa
	rm left.fa right.fa &
	touch $(RUN)/jellyfish.1.finished

$(RUN).Trinity.fasta:$(RUN)_left.fastq $(RUN)_right.fastq
	Trinity --min_kmer_cov $(MINK) --seqType $(SEQ) --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G --bflyCPU $(BCPU) \
	--left $(RUN)_left.fastq --right $(RUN)_right.fastq --group_pairs_distance 999 --CPU $(CPU) --output $(RUN) >>$(RUN).trinity.pe.log





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
