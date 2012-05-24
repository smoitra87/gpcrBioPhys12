"""
Cleans the blast file for Rhodopsin to generate a clean alignment that
we can use
"""

import os, sys
import Bio
from Bio import Seq,SeqRecord,Alphabet,SeqIO,AlignIO

fname = 'blast_94.fasta'
DATA_DIR = 'data'
fpath = os.path.join(DATA_DIR,fname)
outf = os.path.splitext(fname)[0]+'_clean.fasta'
outfpath = os.path.join(DATA_DIR,outf)


seqrs = list(SeqIO.parse(fpath,"fasta"));
for seqr in seqrs : 
	seq = seqr.seq
	while 'X' in seq : 
		pos = seq.find('X')	
		seqm = seq.tomutable()
		seqm[pos] = '-'
		seq = seqm.toseq()
	seqr.seq = seq		


SeqIO.write(seqrs,open(outfpath,'w'),"fasta")
