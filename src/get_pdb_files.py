#! /usr/bin/env python
from Bio.PDB import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import urllib2,sys,os

t = sys.argv[1]
f = open(t,'r')
pdb_list_f = open('../pdb_list/identical_MHC_pdb_list.txt','w')
rec_list = list()
all_pdbs = list()

if not(os.path.isdir('../raw_pdbs/identical_MHC/')):
	os.mkdir('../raw_pdbs/identical_MHC/')

for line in f:
	line = line.strip()
	if line.startswith("queryid"):
		continue
		
	tokens = line.split('\t')
	pdb_id = tokens[1]
	
	pdb_tokens = pdb_id.split('|')
	pdb = pdb_tokens[3]
	chain = pdb_tokens[4]
	
	if os.path.isfile('../raw_pdbs/identical_MHC/'+pdb+'.pdb'):
		print 'file ../raw_pdbs/identital_MHC/'+pdb+'.pdb already exists'

	else:
		print 'getting ../raw_pdbs/identical_MHC/'+pdb+'.pdb' 
		url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdb
  	
  		pdb_file = urllib2.urlopen(url).read()
  		
		out_file = open('../raw_pdbs/identical_MHC/'+pdb+'.pdb', 'w')
		out_file.write(pdb_file)
		out_file.close()
	
	all_pdbs.append(pdb)

unique_pdbs = list(set(all_pdbs))	
parser = PDBParser(QUIET = True)
for pdb in unique_pdbs:
	structure = parser.get_structure(pdb, '../raw_pdbs/identical_MHC/'+pdb+'.pdb')
	ppb=PPBuilder()
	model=structure[0]
	if model.has_id('P'):
		chain = model['P']
		polypeptide = ppb.build_peptides(chain)
		seq = polypeptide[0].get_sequence()
		
		if len(seq) < 12:
			rec = SeqRecord(seq, 
				id = pdb+"_P", 
				description="MHC peptide")
			rec_list.append(rec)
			pdb_list_f.write(pdb+'\n')
		else:
			os.remove('../raw_pdbs/identical_MHC/'+pdb+'.pdb')

	else:
		os.remove('../raw_pdbs/identical_MHC/'+pdb+'.pdb')

pdb_list_f.close()
f.close()

if not(os.path.isdir('../peptide_seqs/')):
	os.mkdir('../peptide_seqs/')

if not(os.path.isdir('../peptide_seqs/unaligned')):
	os.mkdir('../peptide_seqs/unaligned')
	 
SeqIO.write(rec_list, "../peptide_seqs/unaligned/identical_MHC_peptides.fa","fasta")
