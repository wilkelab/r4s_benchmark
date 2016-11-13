import sys

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def main(argv):

	if len(argv) != 3: # wrong number of arguments
		print """Usage:
format_tree_node_id.py <tree_file> <reformated_tree_file>
"""
		sys.exit()

	tree_file=argv[1]
	out_tree_file=argv[2]
	
	f=open(tree_file,"r")
	out=open(out_tree_file,"w")
	for line in f:
		line=line.strip()
		for ch in line:	
		
		
if __name__ == "__main__":
	main(sys.argv)