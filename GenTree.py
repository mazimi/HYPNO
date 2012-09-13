#!/opt/python/bin/python2.7

import os, re

class GenTree:

    # __init__: Constructor method
    # input: self
    # output: none
    def __init__( self ):
        pass

    #Generates tree using FastTree
    def makeTree( self, MSA , name):
    	os.system("FastTree -gtr -nt " + MSA + " > " + name)

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
	    		myRegEx1 = r'(\(' + myLongID1 + r'|' + myLongID1 + r')'
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
		    		while (re.search(myRegEx2, newickSubString) == None) or (findSubtreeString(newickSubString) > 0):
			    		subTreeEnd += 1
			    		#print str(iterSubTree+1) + ', ' + str(subTreeEnd)
			    		#print newickSubString
				    	newickSubString = newickString[subTreeStart:subTreeEnd]
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
	    		newickSubTreeFile = open(ID + '/subtree' + str(iterSubTree+1) + '.ml' )
	    		newickSubTree = newickSubTreeFile.readline()
	    		newickSubTree = newickSubTree.replace(';', '')
	    		newickSubTree = newickSubTree.replace('\n', '')

	    		prunedTree = prunedTree.replace('<' + str(iterSubTree+1) + '>', newickSubTree)
	    	iterSubTree += 1

    	return prunedTree


if __name__ == '__main__':
	sys.exit(1)