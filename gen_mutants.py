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
from BeautifulSoup import BeautifulSoup as SoupClass
import BeautifulSoup

DATA_DIR='data'
refseqf='1F88.fasta'
refseqpath = os.path.join(DATA_DIR,refseqf)
mutf = "mut_oznur.dat"
mutfpath = os.path.join(DATA_DIR,mutf)
mutozf = "mut_oznur.fasta"
mutozfpath = os.path.join(DATA_DIR,mutozf)
mutozfoldf = "mut_oznur.fold"
mutozfpath = os.path.join(DATA_DIR,mutozfoldf)

hgmdhtmlf = 'rhodop_hgmd.html' 
hgmdhtmlfpath = os.path.join(DATA_DIR,hgmdhtmlf)
hgmdseqf = "hgmd.fasta"
hgmdseqfpath = os.path.join(DATA_DIR,hgmdseqf)

aa_dict = {
'Ala':	'A',
'Arg':	'R',
'Asn':	'N',
'Asp':	'D',
'Cys':	'C',
'Glu':	'E',
'Gln':	'Q',
'Gly':	'G',
'His':	'H',
'Ile':	'I',
'Leu':	'L',
'Lys':	'K',
'Met':	'M',
'Phe':	'F',
'Pro':	'P',
'Ser':	'S',
'Thr':	'T',
'Trp':	'W',
'Tyr':	'Y',
'Val':	'V'
}


#----------------------------------------------------------------------
# Oznur data

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

#----------------------------------------------------------------------
# HGMD stuff

with open(hgmdhtmlfpath) as fin : 
	html_doc = fin.read()

soup = SoupClass(html_doc)
rows = soup.findAll('tr')
rows = rows[1:]

mut_list = []

for row in rows :
	elems = list(row.childGenerator())
	mut,pos,ref = elems[2],elems[3],elems[6]
	reflink = ref.find('a')['href']
	mut,pos,ref = map(lambda x: getattr(x,'text'),(mut,pos,ref))		
	pos = int(pos)
	if 'Term' in mut : 
		continue
	mut = aa_dict[mut.split('-')[1]]

	mutseqs = refseqs[:pos-1]+mut+refseqs[pos:]
	mutseqr = SeqRecord(seq=Seq(mutseqs),\
		description="Mutant 1F88 seq",id=str(pos)+mut\
		,name="1F88")
	mut_list.append(mutseqr)
	
SeqIO.write(mut_list,open(hgmdseqfpath,'w'),"fasta")
	





