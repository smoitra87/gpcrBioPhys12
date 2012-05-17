"""
Split the MSA into different phylogenetic categories - 
 * Adrenergic receptors
 * Chemokine receptors
 * Rhodopsin-like receptors

"""

import sys,os

import Bio
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
DO_BLAST = True
BUFFER_SIZE = 30

#----------------------------------------------------------------------
# Intiial steps

if DO_BLAST : 

	# Load the alignment
	aln = AlignIO.read(open(msa_path),"fasta")
	print("Ranga alignment length is %d"%(aln.get_alignment_length()))

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

with open(msa_xml_path,"r") as fin : 
	blast_records = list(NCBIXML.parse(fin))

	
for blast_record in blast_records  :
	pass	




