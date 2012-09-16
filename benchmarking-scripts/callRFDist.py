#!/opt/python/bin/python2.7

from RFDist import RFDist
import sys
import dendropy

#Sample usage: python2.7 callRFDist.py '(A,(B,C))' '((A,B),C)'


def main(tree1Newick, tree2Newick):
	
	RFDistFactory = RFDist()
	print RFDistFactory.calcTreeDist(tree1Newick, tree2Newick)


if __name__ == '__main__':
	tree1Newick = sys.argv[1]
	tree2Newick = sys.argv[2]
	main(tree1Newick, tree2Newick)