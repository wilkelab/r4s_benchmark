import sys

def main(argv):

	if len(argv) != 3: # wrong number of arguments
		print """Usage:
format_tree_id.py <tree_file> <reformated_tree_file>
"""
		sys.exit()

	tree_file=argv[1]
	out_tree_file=argv[2]
	
	f=open(tree_file,"r")
	out=open(out_tree_file,"w")
	
	for line in f:
		line=line.strip()
		if "|" in line:
			line=line.replace("|","_")
		
		new_line=''
		for i in range(len(line)):
			if "." in line[i]:
				if line[i-2]==":":
					new_line += line[i]
				else:
					new_line += "_"
			else:
				new_line += line[i]
		
		out.write(new_line+"\n")
	
if __name__ == "__main__":
	main(sys.argv)