#!/opt/python/bin/python2.7

from time import time
import os
import shutil
import argparse
from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
from kerf import Kerf
import re
import DNA2ProtAlign


#Creates a temporary directory to store data files
def initialize(ID, msa, tree):
	os.mkdir(ID)
	shutil.move(msa, ID)
	shutil.move(tree, ID)

#Calls kerf as a class to split trees into
#subtrees based on given threshold
def makeSubtrees(ID, msa, tree, threshold):
	fileMSA = str(msa)
	fileTree = str(tree)
	threshold = float(threshold)
	outputDir = str(ID)
	kerfTree = Kerf()
	kerfTree.kerfRun( fileTree, fileMSA, threshold, outputDir)
	#OLD: subprocess.call(["python kerf.py --msa-file " + ID + '/' + msa + " --tree-file " + ID + '/' + tree + " --output-directory " + ID + " --threshold " + str(threshold)], shell=True)
	#TODO
	#Should do some error handling
	
#Reads Kerf's msa_map.csv file and extracts the
#number of subtree files that were generated
def getNumSubTrees(ID):
	numSubTrees = 0

	if os.path.isfile(ID + '/msa_map.csv'):
		msaMap = open(ID + '/msa_map.csv', 'rU')
		for line in msaMap:
			match = re.search(r'^(\d+)\,', line)
			if match:
				numSubTrees = int(match.group(1))
				numSubTrees += 1 				#Since Kerf subtree files are indexed from 0
	
	return numSubTrees

#Retreives DNA sequences for each subtree
def getDNASeqs(ID):

	numSubTrees = getNumSubTrees(ID)

	for i in xrange(0, numSubTrees):
		#Retreive Uniprot Accession IDs from subtrees
		#Restricted to Uniprot ID format listed at: http://www.uniprot.org/manual/accession_numbers
		print 'Subtree #' + str(i)	#DEBUG: print subtree number
		listAccessionIDs = ()
		tree = Phylo.read(ID + '/' + str(i) + '.nh', 'newick')
		for clade in tree.find_clades():
			if clade.name:
				match = re.search(r'\|([A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]|[O-Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9])\|', clade.name)
				if match:
					print '\t' + match.group(1)	#DEBUG: print all accession IDs for this subtree
					listAccessionIDs.append(match.group(1))
		#TODO
		#Pass list of accession IDs to DNA sequence lookup script
		

def main():
	#Parse arguments specifying MSA and TREE files
	parser = argparse.ArgumentParser(description='Method for HYbrid Protein NucleOtide phylogenetic tree construction')
	parser.add_argument('--msa-file', required=True)
	parser.add_argument('--tree-file', required=True)
	args = parser.parse_args()
	msa = str(args.msa_file)
	tree = str(args.tree_file)
	threshold = 90.0

	ID = time()									#Used to create directory to store files
	initialize(str(ID), msa, tree)				#Create dir and move input files
	makeSubtrees(str(ID), msa, tree, threshold)	#Run Kerf
	getDNASeqs(str(ID))


if __name__ == '__main__':
	main()