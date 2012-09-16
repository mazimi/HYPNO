#!/usr/bin/env python2.7

#INPUT: Two newick strings
#OUTPUT: RF distance between trees

import os
import dendropy


class RFDist:

	# __init__: Constructor method
	# input: self
	# output: none
	def __init__(self):
		pass

	def calcTreeDist(self, tree1Newick, tree2Newick):

		tree1 = dendropy.Tree.get_from_string(tree1Newick, schema="newick")
		tree2 = dendropy.Tree.get_from_string(tree2Newick, schema="newick")

		return tree1.symmetric_difference(tree2)


if __name__ == '__main__':
	sys.exit(1)