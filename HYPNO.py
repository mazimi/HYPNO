#!/opt/python/bin/python2.7

import os, sys, shutil, argparse, re
from time import time
from Bio import AlignIO, Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from kerf import Kerf
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
	with open(ID+'/HYPNO.debug','w') as debugFh:
		debugFh.write('HYPNO debug info:\n\n')

#Retreives DNA sequences for each subtree
#from Kerf generated csv MSA map file
def getDNASeqs(ID, msa, pred, execMin):
	numSubTrees = 0
	listAccessionIDs = [[]]
	listLongID = [[]]
	treeHierarchy = []

	#Extract Uniprot IDs from Kerf CSV output
	if os.path.isfile(ID + '/' + msa.split('.')[0] + '.csv'):
		msaMap = open(ID + '/' + msa.split('.')[0] + '.csv', 'rU')
		#Retreive Uniprot Accession IDs from subtrees
		#Restricted to Uniprot ID format listed at: http://www.uniprot.org/manual/accession_numbers
		i = 0
		for line in msaMap:
			match = re.search(r'^(\d+).*([A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]|[O-Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9])', line)
			if match:
				if numSubTrees != int(match.group(1)):
					numSubTrees = int(match.group(1))
					listAccessionIDs.append([])
					listLongID.append([])
					treeHierarchy.append(0)
					i += 1
					#print 'Subtree #' + str(numSubTrees) #DEBUG: print subtree number
				#print '\t' + match.group(2)	#DEBUG: print all accession IDs for this subtree
				listAccessionIDs[numSubTrees].append(match.group(2))
				longMatch = re.search(r'^(\d+)\,\s(.*)\n', line)
				listLongID[numSubTrees].append(longMatch.group(2))
				treeHierarchy[i-1] += 1
			else:
				print '** HYPNO input error: given MSA headers must contain valid UniProt accessions.', \
					  '\tAccession missing for given line: '+line
				sys.exit(1)
		#TODO
		#Pass list of accession IDs to DNA sequence lookup script (uniprot.py)
		#Create FASTA files with DNA sequence for each clade (uniprot.py)

	# number of missed / totalTried nucleotide retrievals must be >= execMin
	missed, totalTried, missedList = float(0), float(0), []
	#Retreive DNA sequences using Uniprot IDs (RTFetch)
	for i in xrange(1, numSubTrees+1):
		DNASeqs2Write = []
		nucFetch = RTfetch()
		#print 'Subtree #' + str(i) #DEBUG: print subtree number
		
		j = 0
		for uniprotID in listAccessionIDs[i]:
			with open(ID+'/HYPNO.debug','a') as debugFh:
				totalTried += 1
				nucFetchTuple = nucFetch.getSeqID(uniprotID)
				if nucFetchTuple[4] == 'null':
					missed += 1; missedList.append(uniprotID)
					if nucFetchTuple[10] != 'null':
						print '\t**HYPNO obsolete UniProt entry warning: ', nucFetchTuple[10]
						debugFh.write('\t**HYPNO obsolete UniProt entry warning: '+nucFetchTuple[10])
					else:
						debugFh.write(nucFetchTuple[11])			
					continue
				debugFh.write(nucFetchTuple[11])
				if float(nucFetchTuple[6]) < pred:
					missed += 1; missedList.append(uniprotID)
					debugFh.write('\tPredicted protein sequence for retrieved nucleotide sequence does not match '+ \
									'user provided threshold '+str(pred)+' for percent identity to expected protein '+ \
									'sequence (only '+str(nucFetchTuple[6])+'). Continuing execution.\n')
					continue
				DNASeqs2Write.append(SeqRecord(Seq(nucFetchTuple[9]), id=listLongID[i][j], description="HYPNO FASTA Output"))
				#print '\t' + uniprotID + ' - ' + nucFetchTuple[9] + '\n' #DEBUG: print uniprot ID and sequence
				j += 1

		#Write DNA sequences for each subtree to file
		output_handle = open(ID + '/DNAseqs' + str(i) + '.fasta', 'w')
		SeqIO.write(DNASeqs2Write, output_handle, "fasta")
		output_handle.close()

	percentPassed = 100 * (1 - missed / totalTried)
	if not percentPassed >= execMin:
		print "\n** HYPNO execution error: Nucleotide sequences were retrieved for only "+str(percentPassed)+ \
				" percent of attempted sequences, with the remaining attempts yielding DNA sequences below the "+ \
				str(execMin)+" --n threshold.\nPossible debugging solutions are as follows:\n"+ \
				"\t(a) Lower the --n or --s threshold values to allow more permissive HYPNO execution.\n"+ \
				"\t(b) Rerun the program with a different input tree and MSA, as it may be the case that "+ \
				"no reliable nucleotide sequence exists for a given UniProt accession.\n"+ \
				"Accessions for which nucleotide sequences could not be retrieved were the following: "+str(missedList)+"\n"+ \
				"Please refer to the debug file HYPNO.debug for more information on the queries that were attempted "+ \
				"for each accession."
		sys.exit(1)

	return numSubTrees, treeHierarchy, listLongID

#For all subtrees, takes protein alignment and
#retreived DNA sequence and inserts DNA into
#original protein alignment
#TODO: Determine if DNA should be re-aligned
def alignDNASeqs(ID, msa, tree, numSubTrees):
	for i in xrange(1, numSubTrees+1):
		fileProt = ID + '/' + msa.split('.')[0] + 'sf' + str(i) + '.a2m'
		fileDNA = ID + '/DNAseqs' + str(i) + '.fasta'
		fileOutput = ID + '/' + 'subtree' + str(i) + '.a2m'
		fileAllTree = ID + '/' + tree

		myDNA = DNA2ProtAlign()
		myDNA.alignDNAseqs(fileProt, fileDNA, fileOutput, fileAllTree)

#For all subtrees, takes DNA alignment and
#generates new subtree using GTR algorithm
def makeSubTrees(ID, numSubTrees):
	#TODO: Is A2M = Aligned FASTA? If not, convert. How to handle lower case columns? Delete?
	for i in xrange(1,numSubTrees+1):
		MSA = ID + '/' + 'subtree' + str(i) + '.a2m'
		outputName = ID + '/' + 'subtree' + str(i) + '.ml'
		myTree = GenTree()
		myTree.makeTree(MSA , outputName, ID)


#Takes newick format output of each subtree
#and places back into site where it was
#pruned in original tree
def mergeTree(ID, tree, treeHierarchy, listLongID):
	#TODO: re-calculate lengths (FastTree?)
	tree = ID + '/' + tree

	myTree = GenTree()
	prunedTree = myTree.pruneTree(tree, treeHierarchy, listLongID)
	mergedTree = myTree.insertSubTrees(ID, prunedTree, treeHierarchy)

	output_handle = open(ID + '/hybridTree.ml', "w")
	output_handle.write(mergedTree)
	output_handle.close()

# Checks for proper input formatting
def validateInputs(msa, tree):
	# Check for existence and proper FASTA formatting of input MSA
	try:
		msaHandle = open(msa, "rU")
	except:
		print '** HYPNO input error: Given MSA file location does not exist or is not accessible: '+msa
		sys.exit(1)
	try:
		SeqIO.parse(msaHandle, "fasta").next()
	except:
		print '** HYPNO input error: improper MSA file format, must be aligned FASTA or a2m format: '+msa
		sys.exit(1)	
	try:
		treeHandle = open(tree, "rU")
	except:
		print '** HYPNO input error: Given tree file location does not exist or is not accessible: '+tree
		sys.exit(1)
	try:
		Phylo.read(treeHandle, "newick")
	except:
		print '** HYPNO input error: improper tree file format, must be Newick format: '+msa
		sys.exit(1)
	return 0

def main():
	#Parse arguments specifying MSA and TREE files
	parser = argparse.ArgumentParser(description='Method for HYbrid Protein NucleOtide phylogenetic gene tree reconstruction')
	parser.add_argument('--msa-file', type = str, help='Name of input multiple sequence alignment file (requires aligned FASTA or a2m format)', required=True)
	parser.add_argument('--tree-file', type = str, help='Name of input tree file (requires newick format)', required=True)
	parser.add_argument('--k', default = 90.0, type = float, help='Minimum subtree pairwise percent identity among leaf sequences for subtree topology to be re-estimated using retrieved nucleotide sequences (default: 90.0)', required=False)
	parser.add_argument('--n', default = 95.0, type = float, help='Minimum "predicted protein" percent identity to the expected sequence for a retrieved nucleotide sequence to be accepted (default: 95.0)', required=False)
	parser.add_argument('--s', default = 100.0, type = float, help='Minimum percent of correct retrieved nucleotide sequences (e.g. those passing --n threshold) for program execution to continue (default: 100.0)', required=False)
	args = parser.parse_args()
	msa, tree, kerf, pred, execMin = args.msa_file, args.tree_file, args.k, args.n, args.s 

	ID = time()									#Used to create directory to store files
	print 'Step 0 of 5: Input validation'
	validateInputs(msa,tree)
	initialize(str(ID), msa, tree)				#Create dir and move input files
	print 'Step 1 of 5: Generating subtrees'
	makeSubtrees(str(ID), msa, tree, kerf)		#Run Kerf
	print 'Step 2 of 5: Retrieving DNA sequences'
	numSubTrees, treeHierarchy, listLongID = getDNASeqs(str(ID), msa, pred, execMin)
	print 'Step 3 of 5: Mapping DNA sequences to given protein MSA'
	alignDNASeqs(str(ID), msa, tree, numSubTrees)
	print 'Step 4 of 5: Re-estimating subtree topologies'
	makeSubTrees(str(ID), numSubTrees)
	print 'Step 5 of 5: Reinserting subtrees into gene tree topology'
	mergeTree(str(ID), tree, treeHierarchy, listLongID)
	print 'HYPNO execution completed.'


if __name__ == '__main__':
	main()