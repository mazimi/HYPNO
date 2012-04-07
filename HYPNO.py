#!/opt/python/bin/python2.7

from time import time
import os
import shutil
import argparse
from Bio import AlignIO, Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from kerf import Kerf
import re
from RTfetch import RTfetch
from DNA2ProtAlign import DNA2ProtAlign
from GenTree import GenTree


#Creates a temporary directory to store data files
def initialize(ID, msa, tree):
	os.mkdir(ID)
	shutil.copy(msa, ID)	#TODO: change back to move
	shutil.copy(tree, ID)

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
	listAccessionIDs = [[]]

	#Extract Uniprot IDs from Kerf CSV output
	if os.path.isfile(ID + '/' + msa.split('.')[0] + '.csv'):
		msaMap = open(ID + '/' + msa.split('.')[0] + '.csv', 'rU')
		#Retreive Uniprot Accession IDs from subtrees
		#Restricted to Uniprot ID format listed at: http://www.uniprot.org/manual/accession_numbers
		for line in msaMap:
			match = re.search(r'^(\d+)\,.*\|([A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]|[O-Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9])\|', line)
			if match:
				if numSubTrees != int(match.group(1)):
					numSubTrees = int(match.group(1))
					listAccessionIDs.append([])
					#print 'Subtree #' + str(numSubTrees) #DEBUG: print subtree number
				#print '\t' + match.group(2)	#DEBUG: print all accession IDs for this subtree
				listAccessionIDs[numSubTrees].append(match.group(2))
		#TODO
		#Pass list of accession IDs to DNA sequence lookup script (uniprot.py)
		#Create FASTA files with DNA sequence for each clade (uniprot.py)

	#Retreive DNA sequences using Uniprot IDs (RTFetch)
	for i in xrange(1, numSubTrees+1):
		DNASeqs2Write = []
		nucFetch = RTfetch()
		print 'Subtree #' + str(i) #DEBUG: print subtree number
		
		for uniprotID in listAccessionIDs[i]:
			nucFetchTuple = nucFetch.getSeqID(uniprotID)
			DNASeqs2Write.append(SeqRecord(Seq(nucFetchTuple[9]), id=uniprotID, description="HYPNO FASTA Output"))
			print '\t' + uniprotID + ' - ' + nucFetchTuple[9] + '\n' #DEBUG: print uniprot ID and sequence

		#Write DNA sequences for each subtree to file
		output_handle = open(ID + '/DNAseqs' + str(i) + '.fasta', 'w')
		SeqIO.write(DNASeqs2Write, output_handle, "fasta")
		output_handle.close()

	return numSubTrees

#For all subtrees, takes protein alignment and
#retreived DNA sequence and inserts DNA into
#original protein alignment
#TODO: Determine if DNA should be re-aligned
def alignDNASeqs(ID, msa, numSubTrees):
	for i in xrange(1, numSubTrees+1):
		fileProt = ID + '/' + msa.split('.')[0] + 'sf' + str(i) + '.a2m'
		fileDNA = ID + '/DNAseqs' + str(i) + '.fasta'
		fileOutput = ID + '/' + 'subtree' + str(i) + '.a2m'

		myDNA = DNA2ProtAlign()
		myDNA.alignDNAseqs(fileProt, fileDNA, fileOutput)

#For all subtrees, takes DNA alignment and
#generates new subtree using GTR algorithm
def makeSubTrees(ID, numSubTrees):
	#TODO: Is A2M = Aligned FASTA? If not, convert. How to handle lower case columns? Delete?
	for i in xrange(1,numSubTrees+1):
		MSA = ID + '/' + 'subtree' + str(i) + '.a2m'
		outputName = ID + '/' + 'subtree' + str(i) + '.ml'

		myTree = GenTree()
		myTree.makeTree(MSA , outputName)


#Takes newick format output of each subtree
#and places back into site where it was
#pruned in original tree
def mergeTree(ID, numSubTrees):
	#TODO: will likely need manual re-insertion to do opposite of KERF
	#		re-calculate lengths (FastTree?)

	numSubTrees = 'changethis'


def main():
	#Parse arguments specifying MSA and TREE files
	parser = argparse.ArgumentParser(description='Method for HYbrid Protein NucleOtide phylogenetic tree construction')
	parser.add_argument('--msa-file', required=True)
	parser.add_argument('--tree-file', required=True)
	args = parser.parse_args()
	msa = str(args.msa_file)
	tree = str(args.tree_file)
	threshold = 95.0

	ID = time()									#Used to create directory to store files
	initialize(str(ID), msa, tree)				#Create dir and move input files
	makeSubtrees(str(ID), msa, tree, threshold)	#Run Kerf
	numSubTrees = getDNASeqs(str(ID), msa)
	alignDNASeqs(str(ID), msa, numSubTrees)
	makeSubTrees(str(ID), numSubTrees)
	#mergeTree(str(ID), numSubTrees)


if __name__ == '__main__':
	main()