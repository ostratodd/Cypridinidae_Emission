import sys, getopt
from Bio import Phylo

def main(argv):
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print 'convertTree.py -i <inputfile Newick> -o <outputfile Nexus>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'convertTree.py -i <inputfile> -o <outputfile>'
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
   tree = Phylo.read(inputfile, "newick")
   Phylo.draw_ascii(tree)
   Phylo.write(tree, outputfile, "nexus")
if __name__ == "__main__":
   main(sys.argv[1:])


