#!/opt/python/bin/python2.7

from time import time
import os
import shutil
import argparse
from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
import subprocess
import re


#Creates a temporary directory to store data files
def initialize(ID, msa, tree):
	os.mkdir(ID)
	shutil.move(msa, ID)
	shutil.move(tree, ID)

#Runs kerf.py as a subprocess to split trees into
#subtrees based on given threshold
def makeSubtrees(ID, msa, tree, threshold):
	subprocess.call(["python kerf.py --msa-file " + ID + '/' + msa + " --tree-file " + ID + '/' + tree + " --output-directory " + ID + " --threshold " + str(threshold)], shell=True)
	#Should do some error handling
	#NOTE: Kerf craps out if threshold is above highest PID in tree (ie. no branches to prune, use threshold 100.0 to reproduce). Need to fix this.

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
		#Call script to retreive DNA sequences:
		#ideally we retreive individual DNA sequences
		#rather than that of an entire subtree
		#this can be done in parallel
		print "placeholder"
		

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