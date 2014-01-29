#!/usr/bin/make -rRf

#USAGE:

SHELL=/bin/bash -o pipefail

TRINITY ?= $(shell which 'Trinity.pl')
MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
TRINDIR := $(dir $(firstword $(TRINITY)))
PATH:=$(MAKEDIR):$(PATH)

BCODES=$(MAKEDIR)barcodes.fa
TRIM = 
MEM=10
CPU=6
RUN=run
READ1=left.fastq
READ2=right.fastq
MUS := Mus_musculus.GRCm38.71.cdna.all.fa
PFAM := /media/macmanes/hd2/imitator/pfam/Pfam-AB.hmm.bin
BLAST_DB= $(MAKEDIR)mus_pep
NUM_SUBSAMP=10000000

TRINITY ?= $(shell which 'Trinity.pl')
MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
TRINDIR := $(dir $(firstword $(TRINITY)))
PATH:=$(MAKEDIR):$(PATH)


all:analysis_files/$(RUN).Trinity.fasta analysis_files/$(RUN).Trinity.fasta.pslx analysis_files/$(RUN).Trinity.fasta.pep analysis_files/$(RUN).blast

trin:analysis_files/$(RUN).Trinity.fasta
pslx:analysis_files/$(RUN).Trinity.fasta.pslx
subsamp:raw.$(NUM_SUBSAMP).$(READ1) raw.$(NUM_SUBSAMP).$(READ2)
pep:analysis_files/$(RUN).Trinity.fasta.pep
blast:analysis_files/$(RUN).blast

raw.$(NUM_SUBSAMP).$(READ1) raw.$(NUM_SUBSAMP).$(READ2):$(READ1) $(READ2) 
		@echo TIMESTAMP: `date +'%a %d%b%Y  %H:%M:%S'` ---SUBSAMPLING--- '\n\n'
		python ${MAKEDIR}/subsampler.py $(NUM_SUBSAMP) $(READ1) $(READ2)
		mv subsamp_1.fastq raw.$(NUM_SUBSAMP).$(READ1)
		mv subsamp_2.fastq raw.$(NUM_SUBSAMP).$(READ2)        
analysis_files/$(RUN).Trinity.fasta.pslx:analysis_files/$(RUN).Trinity.fasta
		@echo TIMESTAMP: `date +'%a %d%b%Y  %H:%M:%S'` ---PSLXing--- '\n\n'
		$(TRINDIR)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target ${MAKEDIR}/$(MUS) \
		--query analysis_files/$(RUN).Trinity.fasta; mv *pslx analysis_files/
		rm *summary *maps *selected
analysis_files/$(RUN).Trinity.fasta.pep:analysis_files/$(RUN).Trinity.fasta
		@echo TIMESTAMP: `date +'%a %d%b%Y  %H:%M:%S'` ---TRANSDECODER--- '\n\n'
		$(TRINDIR)/trinity-plugins/transdecoder/TransDecoder --quiet --CPU $(CPU) -t analysis_files/$(RUN).Trinity.fasta \
		--search_pfam $(PFAM); \
		@rm *bed *gff3 *dat *tbl *cds; mv $(RUN).Trinity.fasta.transdecoder.pep analysis_files/$(RUN).Trinity.fasta.pep
		rm -fr transdecoder*
analysis_files/$(RUN).blast:analysis_files/$(RUN).Trinity.fasta
		@echo TIMESTAMP: `date +'%a %d%b%Y  %H:%M:%S'` ---BLASTp--- '\n\n'
		blastp -query analysis_files/$(RUN).Trinity.fasta.pep -db $(BLAST_DB) -evalue 1e-80 -num_threads 8 \
		-outfmt "6 qseqid sacc pident length evalue" > analysis_files/$(RUN).blast

###
###CONDITIONAL EXECUTION
###


analysis_files/$(RUN).Trinity.fasta:raw.$(NUM_SUBSAMP).$(READ2) raw.$(NUM_SUBSAMP).$(READ1)
			if [ ! -d $(RUN) ] ; then mkdir $(RUN) ; fi
			if [ ! -d analysis_files ] ; then mkdir analysis_files ; fi
ifneq ($(TRIM),)
			java -XX:ParallelGCThreads=32 -Xmx$(MEM)g -jar ${MAKEDIR}/trimmomatic-0.32.jar PE \
			-phred33 -threads $(CPU) \
			raw.$(NUM_SUBSAMP).$(READ1) \
			raw.$(NUM_SUBSAMP).$(READ2) \
			$(TRIM).$(NUM_SUBSAMP).PP.$(READ1) \
			$(TRIM).$(NUM_SUBSAMP).UP.$(READ1) \
			$(TRIM).$(NUM_SUBSAMP).PP.$(READ2) \
			$(TRIM).$(NUM_SUBSAMP).UP.$(READ2) \
			ILLUMINACLIP:$(BCODES):2:40:15 \
			LEADING:$(TRIM) \
			TRAILING:$(TRIM) \
			SLIDINGWINDOW:4:$(TRIM) \
			MINLEN:25 2>> trim10.log;
			cat $(TRIM).$(NUM_SUBSAMP).PP.$(READ1) $(TRIM).$(NUM_SUBSAMP).UP.$(READ1) $(TRIM).$(NUM_SUBSAMP).UP.$(READ2) > $(RUN).left.$(TRIM).fq ;
			mv $(TRIM).$(NUM_SUBSAMP).PP.$(READ2) $(RUN).right.$(TRIM).fq;
#trim
			$(TRINDIR)/Trinity.pl --bflyGCThread 25 --full_cleanup --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G \
			--left $(RUN).left.$(TRIM).fq --right $(RUN).right.$(TRIM).fq --group_pairs_distance 999 \
			--CPU $(CPU) --CuffFly --output $(RUN);
			mv $(RUN).Trinity.fasta analysis_files/
else

#fastool
			@echo TIMESTAMP: `date +'%a %d%b%Y  %H:%M:%S'` ---FASTOOL--- \n\n
			fastool --illumina-trinity --to-fasta raw.$(NUM_SUBSAMP).$(READ1) > $(RUN)/left.fa 2> $(RUN)/reads.left.fq.readcount
			fastool --illumina-trinity --to-fasta raw.$(NUM_SUBSAMP).$(READ2) > $(RUN)/right.fa 2> $(RUN)/reads.right.fq.readcount
			cat $(RUN)/left.fa $(RUN)/right.fa > $(RUN)/both.fa
			touch $(RUN)/both.fa.read_count
#jellyfish
			@echo TIMESTAMP: `date +'%a %d%b%Y  %H:%M:%S'` ---JELLYFISH--- \n\n
			~/jellyfish-2.0.0/bin/jellyfish count -m25 -t$(CPU) -C -s2G -o $(RUN)/jf.5prime <(python ${MAKEDIR}/5prime.py $(RUN)/both.fa)
			~/jellyfish-2.0.0/bin/jellyfish bc -m 25 -s3G -C -t $(CPU) -o  $(RUN)/both.bf  $(RUN)/both.fa
			~/jellyfish-2.0.0/bin/jellyfish count -m25 -t$(CPU) -C -s2G --bc $(RUN)/both.bf -o $(RUN)/jf.3prime <(python ${MAKEDIR}/3prime.py $(RUN)/both.fa)
			~/jellyfish-2.0.0/bin/jellyfish  merge $(RUN)/jf.3prime $(RUN)/jf.5prime -o $(RUN)/merged.jf 
			~/jellyfish-2.0.0/bin/jellyfish dump -o $(RUN)/jellyfish.kmers.fa $(RUN)/merged.jf
			touch $(RUN)/jellyfish.1.finished
#trin	
			$(TRINDIR)/Trinity.pl --bflyGCThread 25 --full_cleanup --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G \
			--left raw.$(NUM_SUBSAMP).$(READ1) --right raw.$(NUM_SUBSAMP).$(READ2) --group_pairs_distance 999 \
			--CPU $(CPU) --CuffFly --output $(RUN);
			mv $(RUN).Trinity.fasta analysis_files/
endif
