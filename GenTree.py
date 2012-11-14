#!/opt/python/bin/python2.7

import re, subprocess
from subprocess import Popen
from Bio import Phylo

class GenTree:

	# __init__: Constructor method
	# input: self
	# output: none
	def __init__( self ):
		pass

	#Generates tree using FastTree and performs midpoint rooting
	def makeTree( self, MSA , name, ID):
		cmd = "FastTree -gtr -nt " + MSA + " > " + name
		pro = Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

		#Perform midpoint root
		#tree = Phylo.read(name, 'newick')
		#tree.root_at_midpoint()
		#Phylo.write(tree, name, 'newick')

		with open(ID+'/HYPNO.debug','a') as debugFh:
			for line in pro.communicate():
				debugFh.write(line.rstrip())
			debugFh.write('\n')
		return 0

	#Recalculate branch lengths keeping topology fixed
	def makeTreeBranchLengths( self, tree, MSA , name, ID):
		#Filter taxa in MSA that are not present in final tree
		#TODO track number of accessions in MSA and Tree and report mismatch
		def filterMSA( tree, MSA, ID):
			msaFile = open(MSA, 'rU')
			treeFile = open(tree, 'rU')
			msaFilteredFile = open(MSA + '_filtered', 'w')
			flagKeepSeq = False

			treeString = treeFile.read()

			for line in msaFile:
				matchHeader = re.search(r'^\>', line)
				if matchHeader:
					matchUniprot = re.search(r'([A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]|[O-Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9])', line)
					if matchUniprot:
						matchUniprotTree = re.search(re.escape(matchUniprot.group(1)), treeString)
						if matchUniprotTree:
							flagKeepSeq = True
						else:
							flagKeepSeq = False
					else:
						flagKeepSeq = False

				if flagKeepSeq:
					msaFilteredFile.write(line)

			msaFile.close()
			treeFile.close()
			msaFilteredFile.close()
			return 0

		filterMSA(tree, MSA, ID)
		cmd = "FastTree -nome -mllen -intree " + tree + " " + MSA + "_filtered > " + name
		pro = Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
		with open(ID+'/HYPNO.debug','a') as debugFh:
			for line in pro.communicate():
				debugFh.write(line.rstrip())
			debugFh.write('\n')
		return 0


	#Takes in a tree and prunes the larger tree where there exist
	#identified subtrees with 3 or more leafs (based on Kerf CSV).
	def pruneTree( self, tree, treeHierarchy, listLongIDs):

		#Check to see that the number of open paranthesis
		#and close paranthesis are balanced
		#(returns the difference between the number of
		#open and close paranthesis in a string)
		def findSubtreeString(newickSubString):
			countOpen = newickSubString.count('(')
			countClose = newickSubString.count(')')
			return (countOpen - countClose)

		newickHandle = open(tree, 'rU')
		newickString = newickHandle.readline()

		iterSubTree = 0
		for subtree in treeHierarchy:
			if subtree >= 3:
				#Build dynamic regex string
				myLongID1 = listLongIDs[iterSubTree+1][0]
				myLongID1 = myLongID1.replace('|', '\|')
				#myRegEx1 = r'(\(' + myLongID1 + r'|' + myLongID1 + r')'
				myRegEx1 = r'(' + myLongID1 + r')'
				#Check for first accession ID within substring
				match1 = re.search(myRegEx1, newickString)
				if match1:
					#Pick out substring from match to end of newick string
					subTreeStart = match1.start()
					subTreeEnd = match1.end()
					newickSubString = newickString[subTreeStart:subTreeEnd]
					#Search for last accession ID in substring starting from
					#first match and expanding string length until number of
					#paranthesis is balanced and last accession is found
					myLongID2 = listLongIDs[iterSubTree+1][treeHierarchy[iterSubTree]-1]
					myLongID2 = myLongID2.replace('|', '\|')
					myRegEx2 = r'(' + myLongID2 + r')'
					#print myRegEx1 + ', ' + myRegEx2
					while (re.search(myRegEx2, newickSubString) == None): #or (findSubtreeString(newickSubString) > 0):
						subTreeEnd += 1
						#print str(iterSubTree+1) + ', ' + str(subTreeEnd)
						#print newickSubString
						newickSubString = newickString[subTreeStart:subTreeEnd]
					print newickSubString
					#Check for nesting in parentheses
					balanceParentheses = findSubtreeString(newickSubString)
					while (balanceParentheses != 0):
						if (balanceParentheses > 0):
							subTreeEnd += 1
						elif (balanceParentheses < 0):
							subTreeStart -= 1
						newickSubString = newickString[subTreeStart:subTreeEnd]
						balanceParentheses = findSubtreeString(newickSubString)
						print subTreeStart, subTreeEnd, newickSubString
					#Grab additional parentheses belonging to this clade
					subTreeStart -= 1
					subTreeEnd += 1
					newickSubString = newickString[subTreeStart:subTreeEnd]
					while (newickSubString[0] == '(' ) and (newickSubString[-1] == ')' ):
						subTreeStart -= 1
						subTreeEnd += 1
						newickSubString = newickString[subTreeStart:subTreeEnd]
					#Strip out non-parentheses
					newickSubString = newickString[subTreeStart+1:subTreeEnd-1]
					print subTreeStart, subTreeEnd, newickSubString
					#Once subtree is found, replace it with the subtree number placed in brackets
					newickString = newickString.replace(newickSubString, '<' + str(iterSubTree+1) + '>')
			iterSubTree += 1
		
		return newickString

	#Goes through subtree and finds locations of subtrees (number within bracket)
	#and replaces them with the newly generated subtrees from FastTree
	def insertSubTrees( self, ID, prunedTree, treeHierarchy):
		iterSubTree = 0
		for subtree in treeHierarchy:
			if subtree >= 3:
				#Perform midpoint root
				tree = Phylo.read(ID + '/subtree' + str(iterSubTree+1) + '.ml', 'newick')
				tree.root_at_midpoint()
				Phylo.write(tree, ID + '/subtree' + str(iterSubTree+1) + '.ml', 'newick')

				newickSubTreeFile = open(ID + '/subtree' + str(iterSubTree+1) + '.ml' )
				newickSubTree = newickSubTreeFile.readline()
				newickSubTree = newickSubTree.replace(';', '')
				newickSubTree = newickSubTree.replace('\n', '')

				prunedTree = prunedTree.replace('<' + str(iterSubTree+1) + '>', newickSubTree)
			iterSubTree += 1

		return prunedTree


if __name__ == '__main__':
	sys.exit(1)