#! /usr/bin/env python
###This script takes in a virus name and search for all available protein structure in protein data bank that match 'HLA-A2 <virus_name>'.
from Bio.PDB import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import urllib2,sys,os

virus = sys.argv[1] 
list_file = virus+'_pdb_list.txt'

if not(os.path.isdir('../pdb_list/')):
	os.mkdir('../pdb_list/')
pdb_list_f = open('../pdb_list/'+list_file,'w')

url = 'http://www.rcsb.org/pdb/rest/search'

queryText = """

<?xml version="1.0" encoding="UTF-8"?>
<orgPdbQuery>

<queryType>org.pdb.query.simple.AdvancedKeywordQuery</queryType>

<description>Text Search for: HLA-A2 %s</description>

<keywords>HLA-A2 %s</keywords>

</orgPdbQuery>

""" %(virus,virus)

print "query:\n", queryText

print "querying PDB...\n"

req = urllib2.Request(url, data=queryText)

f = urllib2.urlopen(req)

result = f.read()

if result:

    print "Found number of PDB entries:", result.count('\n')
    result = result.rstrip()
    raw_pdb_id_list = result.split('\n')
    print "PDB entries found: ", "  ".join(raw_pdb_id_list)
    
else:

    print "Failed to retrieve results"

if not(os.path.isdir('../raw_pdbs/'+virus)):
	os.mkdir('../raw_pdbs/'+virus) 

print "Downloading PDB files found ... "

filtered_pdb_id_list = list()
for pdb in raw_pdb_id_list:  
	url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdb
  	
  	pdb_file = urllib2.urlopen(url).read()
  	title_count=pdb_file.find("TITLE")
  	pdb_title=pdb_file[title_count:title_count+200]
  	
  	if virus in pdb_title:
  		pdb_list_f.write(pdb+'\n')
  		filtered_pdb_id_list.append(pdb)
  		if os.path.isfile('../raw_pdbs/'+virus+'/'+pdb+'.pdb'):
  			print 'file ../raw_pdbs/'+virus+'/'+pdb+'.pdb already exists'
			continue
		else:
			print 'writing ../raw_pdbs/'+virus+'/'+pdb+'.pdb' 
			out_file = open('../raw_pdbs/'+virus+'/'+pdb+'.pdb', 'w')
			out_file.write(pdb_file)
			out_file.close()

rec_list = list()
for pdb in filtered_pdb_id_list:
	parser = PDBParser(QUIET = True)

	structure = parser.get_structure(pdb, '../raw_pdbs/'+virus+'/'+pdb+'.pdb')
	ppb=PPBuilder()
	polypeptide = ppb.build_peptides(structure[0]['P'])
	seq = polypeptide[0].get_sequence()
	rec = SeqRecord(seq, 
		id = pdb+"_P", 
		description=virus+" peptide")
	rec_list.append(rec) 

if not(os.path.isdir('../peptide_seqs/')):
	os.mkdir('../peptide_seqs/')

if not(os.path.isdir('../peptide_seqs/unaligned')):
	os.mkdir('../peptide_seqs/unaligned')
	 
SeqIO.write(rec_list, "../peptide_seqs/unaligned/"+virus+"_peptide.fa","fasta")
