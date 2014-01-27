#!/usr/bin/python
from Bio import SeqIO
import sys
import os


for seq in SeqIO.parse(open(sys.argv[1], "rU"), "fasta"):
	print "%s%s%s%s" % ('>',seq.id, '\n', seq.seq[-25:])
