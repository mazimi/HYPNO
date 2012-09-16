#!/usr/bin/env python2.7

#INPUT: Newick string (without branch lengths)
#OUTPUT: Generates all possible newick strings in which trifurcations are forced
#		 to be bifurcations
#OUTPUT FORMAT: Tuple (number of new newick strings, list of new newick strings)

import os
import re
import dendropy


class Multifurcations:

	# __init__: Constructor method
	# input: self
	# output: none
	def __init__(self):
		pass

	def genAllBifurcations(self, myNewick):

		#Builds newick string with square brackets around upper most node level branches
		def getMultifurcatingNode(myNode, isRoot):
			#If node is internal, recursively retreive child nodes and build up newick string
			if myNode.is_internal() and isRoot==False:
				myNewickNode = '('
				for childNode in myNode.child_nodes():
					myNewickNode += getMultifurcatingNode(childNode, False) + ','
				myNewickNode = myNewickNode.rstrip(',')
				myNewickNode += ')'
				return myNewickNode
			#In the case that the node is the root node that was submitted, place square
			# brackets around three groups for easy parsing.
			elif myNode.is_internal() and isRoot==True:
				myNewickNode = '('
				for childNode in myNode.child_nodes():
					myNewickNode += '[' + getMultifurcatingNode(childNode, False) + '],'
				myNewickNode = myNewickNode.rstrip(',')
				myNewickNode += ')'
				return myNewickNode
			#If nodes are leafs then return labels
			else:
				return myNode.as_newick_string()

		#Generates three possible combinations of bifurcating trees from a single trifurcating tree
		def expandMultifurcating(myNewickNode):
			bifurcatingNewick = []
			#Square brackets surround three branches at the upper most node
			match = re.search(r'^\(\[(\S+)\],\[(\S+)\],\[(\S+)\]\)$', myNewickNode)

			bifurcatingNewick.append('(' + match.group(1) + ',(' + match.group(2) + ',' + match.group(3) + '))')
			bifurcatingNewick.append('(' + match.group(2) + ',(' + match.group(1) + ',' + match.group(3) + '))')
			bifurcatingNewick.append('(' + match.group(3) + ',(' + match.group(1) + ',' + match.group(2) + '))')

			return bifurcatingNewick



		treeList = []
		treeListToRemove = []
		treeListNewick = []

		treeOrig = dendropy.Tree.get_from_string(myNewick, schema="newick")
		treeList.append(treeOrig)

		#Count number of trees that are expected
		numTreesFinal = 1
		for node in treeOrig.postorder_node_iter():
			if len(node.child_nodes()) > 2:
				if len(node.child_nodes()) > 3:
					#Can't handle more than trifurcation at the moment
					return (0,['I can only handle trifurcating nodes, y\'all.'])
				else:
					#If trifurcation multiply number of final trees by three
					numTreesFinal *= 3

		#Start from children node and traverse up tree looking for trifurcations
		# If a trifurcation is found, original tree is deleted and three new trees
		# are added in which trifurcation is replaced with bifurcations. Process is
		# repeated until all trifurcations are replaced with bifurcations.
		for tree in treeList:
			for node in tree.postorder_node_iter():
				if len(node.child_nodes()) > 2:
					multifurcatingNewick = getMultifurcatingNode(node, True)
					multifurcatingNewick = multifurcatingNewick.replace('[', '')
					multifurcatingNewick = multifurcatingNewick.replace(']', '')
					listBifurcating = expandMultifurcating(getMultifurcatingNode(node, True))
					for newNode in listBifurcating:
						newTree = tree.as_newick_string().replace(multifurcatingNewick, newNode)
						treeList.append(dendropy.Tree.get_from_string(newTree, schema="newick"))
					treeListToRemove.append(tree)
					break

		#This is the step where trees containing trifurcations are removed
		for tree in treeListToRemove:
			treeList.remove(tree)

		#Convert tree objects to newick strings
		for tree in treeList:
			treeListNewick.append(tree.as_newick_string())

		#Make sure that the number of generated trees matches number expected
		if len(treeListNewick) == numTreesFinal:
			return (numTreesFinal, treeListNewick)
		else:
			return (0, ['Number of trees returned doesn\'t match expected.'])

if __name__ == '__main__':
	sys.exit(1)