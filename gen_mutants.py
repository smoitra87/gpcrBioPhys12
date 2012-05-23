""" 

Generates mutants of GPCR sequences and finally outputs an alignment
1U19 is used as the reference sequence

"""

import os,sys,re

import Bio
from Bio import SeqIO,Alphabet,AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from operator import itemgetter


DATA_DIR='data'
refseqf='1F88.fasta'
refseqpath = os.path.join(DATA_DIR,refseqf)
mutf = "mut_oznur.dat"
mutfpath = os.path.join(DATA_DIR,mutf)
mutozf = "mut_oznur.fasta"
mutozfpath = os.path.join(DATA_DIR,mutozf)
mutozfoldf = "mut_oznur.fold"
mutozfpath = os.path.join(DATA_DIR,mutozfoldf)

# Load the refernce sequence
refseqr = list(SeqIO.parse(refseqpath,"fasta"))[0]
refseqs = refseqr.seq.tostring()

# Get the mutant list
with open(mutfpath) as fin : lines = fin.readlines()

lines =  map(str.split,map(str.strip,lines))

for line in lines : 
	mut = line[0]
	line[0] = re.match('(\D)(\d+)(\D)',mut).groups()

lines = filter(lambda l : l[0][2] != 'C',lines)

print("Total number of Oznur mutant sequences: %d"%(len(lines)))

mut_list = []

for line in lines : 
 	prev,pos,mut = line[0]
	pos = int(pos)
	mutseqs = refseqs[:pos-1]+mut+refseqs[pos:]
	mutseqr = SeqRecord(seq=Seq(mutseqs),\
		description="Mutant 1F88 seq",id="".join(line[0])\
		,name="1F88")
	mut_list.append(mutseqr)
	
SeqIO.write(mut_list,open(mutozfpath,'w'),"fasta")


with open(mutozfpath,'w') as fout : 
	fout.write("\n".join(map(itemgetter(1),lines))+"\n")




