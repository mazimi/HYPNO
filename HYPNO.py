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

#Retreives DNA sequences for each subtree
#from Kerf generated csv MSA map file
def getDNASeqs(ID, msa):
	numSubTrees = 0
	listAccessionIDs = []

	if os.path.isfile(ID + '/' + msa.split('.')[0] + '.csv'):
		msaMap = open(ID + '/' + msa.split('.')[0] + '.csv', 'rU')
		#Retreive Uniprot Accession IDs from subtrees
		#Restricted to Uniprot ID format listed at: http://www.uniprot.org/manual/accession_numbers
		for line in msaMap:
			match = re.search(r'^(\d+)\,.*\|([A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]|[O-Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9])\|', line)
			if match:
				if numSubTrees != int(match.group(1)):
					numSubTrees = int(match.group(1))
					print 'Subtree #' + str(numSubTrees) #DEBUG: print subtree number
				print '\t' + match.group(2)	#DEBUG: print all accession IDs for this subtree
				listAccessionIDs.append(match.group(2))
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
	getDNASeqs(str(ID), msa)


if __name__ == '__main__':
	main()