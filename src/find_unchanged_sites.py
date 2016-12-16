import sys
from Bio import AlignIO

def find_unchanged_sites(aln,out):
	total_col=len(aln[0])
	
	for i in range(total_col):
		col = aln[:,i]
		if col == len(col) * col[0]:
			line=str(i+1)+'\tTRUE\n' 
		else:
			line=str(i+1)+'\tFALSE\n' 
		out.write(line)

def main(argv):

	if len(argv) != 3: # wrong number of arguments
		print """Usage:
python find_unchanged_sites.py <aa_aln_file> <output_file>
"""
		sys.exit()

	aln_file = sys.argv[1] 
	out_table = sys.argv[2]

	aln = AlignIO.read(aln_file, "fasta") 
	
	out=open(out_table,"w")
	out.write("pos\tunchanged_site\n")
	find_unchanged_sites(aln,out)

if __name__ == "__main__":
	main(sys.argv)