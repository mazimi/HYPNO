#!/opt/python/bin/python2.7

from Multifurcations import Multifurcations
import sys
import dendropy

#Sample usage: python2.7 callMultifurcations.py '(A,B,C)'


def main(myNewick):
	
	#Following two lines are all you need to use factory, remaining lines are for visualization
	multifurcationFactory = Multifurcations()
	(numTrees, treeListNewick) = multifurcationFactory.genAllBifurcations(myNewick)
	#Above line passes in a string 'myNewick' to the factory and gets back a tuple:
	# numTrees: number of bifurcating trees generated
	#	if numTrees==0, there was an error. Error message is written to treeListNewick
	# treeListNewick: list of all the new bifurcating trees in newick format

	treeOrig = dendropy.Tree.get_from_string(myNewick, schema="newick")
	print 'Original tree: ' + treeOrig.as_newick_string()
	treeOrig.print_plot()

	print 'Generated ' + str(numTrees) + ' bifurcating trees:'
	for treeNewick in treeListNewick:
		tree = dendropy.Tree.get_from_string(treeNewick, schema="newick")
		tree.print_newick()
		tree.print_plot()

if __name__ == '__main__':
	myNewick = sys.argv[1]
	main(myNewick)