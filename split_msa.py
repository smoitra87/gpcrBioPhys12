"""
Split the MSA into different phylogenetic categories - 
 * Adrenergic receptors
 * Chemokine receptors
 * Rhodopsin-like receptors

WARNING :  If BLASTING again for some reason then, run the blast once 
and without blast again. There is some issue with alignment being 
overwritten during blasting. 


"""

import sys,os

import numpy as np

import Bio
from Bio import Alphabet
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import AlignInfo
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

#----------------------------------------------------------------------
# Definitions

DATA_DIR = "data/"
msa_file = "GPCR_ranga.fasta" 
msa_path = os.path.join(DATA_DIR,msa_file)
if not os.path.exists(msa_path) : 
	raise RuntimeError("GPCR ranga file does not exist")
msa_xml_file = "GPCR_ranga_blast.xml"
msa_xml_path = os.path.join(DATA_DIR,msa_xml_file)
E_VALUE = 10e-3
HITLIST_SIZE = 1
DO_BLAST = False
DO_FILTER = False
BUFFER_SIZE = 30
adren_words = ('adrenergic','adrenoceptor')
chemo_words = ('chemokine','cytokine')

#----------------------------------------------------------------------
# Intiial steps

# Load the alignment
aln = AlignIO.read(open(msa_path),"fasta",\
	alphabet=Alphabet.IUPAC.protein)
print("Ranga alignment length is %d"%(aln.get_alignment_length()))


if DO_BLAST : 

	for buffer_id in range(0,len(aln),BUFFER_SIZE):
		aln_buff = aln[buffer_id:(buffer_id + BUFFER_SIZE)]
	
		fas_seq = ""
		for seqr in aln_buff : 
			seqr.seq = seqr.seq.strip('.')
			fas_seq = fas_seq + seqr.format("fasta")
		
		print("Trying to blast Many swquences")
		res_handle = NCBIWWW.qblast(program="blastp",database="nr",\
			sequence=fas_seq,hitlist_size=HITLIST_SIZE,\
			expect=E_VALUE)
		
		# Save results to XML file 
		msa_fname = os.path.splitext(msa_xml_file)[0]+"_"+\
			str(buffer_id)+".xml"
		msa_xml_fpath = os.path.join(DATA_DIR,msa_fname)
		print(msa_xml_fpath)

		with open(msa_xml_fpath,"w") as fout :
			fout.write(res_handle.read())
		res_handle.close()


#----------------------------------------------------------------------
# Parse XML file and classify sequences by looking at top hit def

blast_records = []

for buffer_id in range(0,len(aln),BUFFER_SIZE) :

	# Get name of buffer file
	msa_fname = os.path.splitext(msa_xml_file)[0]+"_"+\
		str(buffer_id)+".xml"
	msa_xml_fpath = os.path.join(DATA_DIR,msa_fname)
	
	with open(msa_xml_fpath,"r") as fin : 
		blast_buff_records = list(NCBIXML.parse(fin))
	
	blast_records.extend(blast_buff_records)
	
# Check all records have been blasted
assert len(blast_records) == len(aln)

seq_class = {}


for blast_record in blast_records  :
	# Check that the blast successfully completed
	assert len(blast_record.alignments) > 0

	# Get the description of blast match sequence
	blast_aln = blast_record.alignments[0]
	hit_def = blast_aln.hit_def
	
	# query id 
	qid = blast_record.query

	#seq_class has seq to class mapping
	
	# check for keywords
 	if reduce(lambda x,y : x or y,\
		map(lambda(w): w in hit_def,adren_words)) :
		seq_class[qid] = "adrenergic"
	elif reduce(lambda x,y : x or y,\
		map(lambda(w): w in hit_def,chemo_words)) :
		seq_class[qid] = "chemokine"
	else : 
		seq_class[qid] = "rhodopsin"
	
# Initialize empty lists for storing sequence records
adr_msa,chemo_msa,rhod_msa = [],[],[]

# Write out MSA
for seqr in aln : 
	if seq_class[seqr.name] == "adrenergic" : 
		adr_msa.append(seqr)
	elif seq_class[seqr.name] == "chemokine" :
		chemo_msa.append(seqr)
	elif seq_class[seqr.name] == "rhodopsin" : 
		rhod_msa.append(seqr)
	else : 
		raise ValueError("seq record name unknown %s"%(seqr.name))

def calc_entropy(counts) : 
	""" Calculates entropy based on a vector of counts"""
	sum_c = sum(counts)


#----------------------------------------------------------------------
# Filter columns if necessary

if DO_FILTER :
	
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
	#	assert False
	#	if e > ENTROPY_THRESH : 
	#		info_cols.append(col_id)
		info_cols.append(e)

msa_path = os.path.join(DATA_DIR,"adrn_msa.fasta")
SeqIO.write(adr_msa,msa_path,"fasta")

msa_path = os.path.join(DATA_DIR,"chemo_msa.fasta")
SeqIO.write(chemo_msa,msa_path,"fasta")

msa_path = os.path.join(DATA_DIR,"rhod_msa.fasta")
SeqIO.write(rhod_msa,msa_path,"fasta")


