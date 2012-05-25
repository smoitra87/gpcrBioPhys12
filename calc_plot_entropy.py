"""
Loads a list of GPCR alignments. Calculates the entropy of the 
alignments. Creates a hist plot of the entropy scores and saves the image

"""
import os,sys,glob

import numpy as np
import pylab as pl

import Bio
from Bio import Seq,SeqIO,AlignIO,SeqRecord,Alphabet

DATA_DIR='data'
#list_msa = glob.glob(os.path.join(DATA_DIR,'*.fasta'))
list_msa = ['GPCR_ranga.fasta','chemo_msa.fasta','adrn_msa.fasta']
list_msa = map(lambda x:os.path.join(DATA_DIR,x),list_msa)

for msa_path in list_msa :
	aln = AlignIO.read(open(msa_path),"fasta",\
		alphabet=Alphabet.IUPAC.protein)
	print("Alignment name is %s"%msa_path)
	print("alignment length is %d"%(aln.get_alignment_length()))

	# Calculate Entropies
	info_cols = [] # Non conserved columns
	# Filter columns according to entropy
	for col_id in range(aln.get_alignment_length()) :
		col = aln.get_column(col_id) 
		ftab = {} # Freq table
		for c in col :
			ftab[c] = ftab.get(c,0) + 1
		counts = ftab.values()
		sum_c = sum(counts)
		p = map(lambda x : (0.0+x)/sum_c,counts) # probabilities
		e = -reduce(lambda x,y:x+y,map(lambda t:t[0]*t[1],\
			zip(p,np.log(p))))
		info_cols.append(e)
	pl.figure()
	pl.hist(info_cols,bins=30)
	pl.title(os.path.splitext(msa_path)[0])
	pl.savefig(os.path.splitext(msa_path)[0]+'_ent.png')
	pl.close()

