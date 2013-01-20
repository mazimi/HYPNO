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
	def makeTreeBranchLengths( self, tree, MSA , name, ID, MSAtype):
		#Filter taxa in MSA that are not present in final tree
		#TODO track number of accessions in MSA and Tree and report mismatch
		def filterMSA( tree, MSA, ID):
			msaFile = open(MSA, 'rU')
			treeFile = open(tree, 'rU')
			msaFilteredFile = open(MSA + '_filtered', 'w')
			flagHeader = False
			flagKeepSeq = False

			treeString = treeFile.read()

			for line in msaFile:
				matchHeader = re.search(r'^\>', line)
				if matchHeader:
					flagHeader = True
					matchUniprot = re.search(r'([A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]|[O-Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9])', line)
					if matchUniprot:
						matchUniprotTree = re.search(re.escape(matchUniprot.group(1)), treeString)
						if matchUniprotTree:
							flagKeepSeq = True
						else:
							flagKeepSeq = False
					else:
						flagKeepSeq = False
				else:
					flagHeader = False

				if (flagHeader == False) and flagKeepSeq:
					msaFilteredFile.write(line)
				elif flagHeader and flagKeepSeq:
					msaFilteredFile.write('>' + matchUniprot.group(1) + ' HYPNO filtered output\n')

			msaFile.close()
			treeFile.close()
			msaFilteredFile.close()
			return 0

		filterMSA(tree, MSA, ID)
		if (MSAtype == 'nuc'):
			cmd = "FastTree -nome -mllen -gtr -nt -intree " + tree + " " + MSA + "_filtered > " + name
		else:
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
		#Read in newick tree
		objTree = Phylo.read(tree, "newick")
		newickString = objTree.format("newick")

		iterSubTree = 0
		for subtree in treeHierarchy:
			if subtree >= 3:
				#Parse out first and last leaves in tree
				myLongID1 = listLongIDs[iterSubTree+1][0]
				myLongID2 = listLongIDs[iterSubTree+1][treeHierarchy[iterSubTree]-1]
				
				#Identify clade containing most recent common ancestor of two leaves
				#This can be modified to contain more than 2 leaves
				objClade = objTree.common_ancestor({"name": myLongID1}, {"name": myLongID2})
				#Convert clade object to tree object
				objSubtree = objTree.from_clade(objClade)
				newickSubString = objSubtree.format("newick")
				newickSubString = newickSubString.rstrip(';\n')
				#Replace substring in 
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