import sys

def main(argv):

	if len(argv) != 3: # wrong number of arguments
		print """Usage:
format_tree_id.py <tree_file> <reformated_tree_file>
"""
		sys.exit()

	fasta_file=argv[1]
	out_fasta_file=argv[2]
	
	f=open(fasta_file,"r")
	out=open(out_fasta_file,"w")
	for line in f:
		if "|" in line:
			line=line.replace("|","_")
		out.write(line)
	
if __name__ == "__main__":
	main(sys.argv)