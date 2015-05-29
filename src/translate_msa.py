##This script takes in nucleotide msa (in fasta format) and converts it to amino acid msa (also in fasta format)
import sys
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

aln_file = sys.argv[1] 
nuc_aln = AlignIO.read(aln_file, "fasta") 

aa_aln = []
for s in nuc_aln:
	aa_seq = s.seq.translate() #translate nucleotide seq to amino acids, returns a Seq object
	aa_seq_r = SeqRecord(aa_seq) #create a SeqRecord object from amino acid Seq object to start a list 
	aa_seq_r.id = s.id #set SeqRecord.id to be printed in the new fasta file
	aa_seq_r.description = '' #set SeqRecord.description to an empty string to avoid <unknown description> being printed in the new fasta file 
	aa_aln.append(aa_seq_r) 

msa = MultipleSeqAlignment(aa_aln) # create an MSA object
base_name = aln_file.split("/")[-1]
AlignIO.write(msa, "../test_msa/aa_msa/"+base_name, "fasta") 