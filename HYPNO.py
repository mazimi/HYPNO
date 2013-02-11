#!/usr/bin/env python2.7

import os, sys, shutil, argparse, re, urllib2
from time import time
from Bio import AlignIO, Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from kerf import Kerf
from RTfetch import RTfetch
from DNA2ProtAlign import DNA2ProtAlign
from GenTree import GenTree

#Creates a temporary directory to store data files
def initialize(ID, msa, tree=None):
	os.mkdir(ID)
	leafLabels, seqHeaders = sanitizeInputs(ID, msa, tree)
	with open(ID+'/HYPNO.debug','w') as debugFh:
		debugFh.write('HYPNO debug info:\n\n')
	return leafLabels, seqHeaders

#Ensure that all entries contain valid Uniprot Accessions
#Produce new input files containing only Uniprot Accessions as identifiers in 'ID'
def sanitizeInputs(ID, msa, tree=None):
	#Sanitize MSA
	leafLabels, seqHeaders = {}, {}
	with open(msa, 'rU') as msaFile:
		with open(ID + '/' + msa, 'w') as msaSanitized:
			lineCount = 0

			for line in msaFile:
				lineCount += 1
				matchHeader = re.search(r'^\>(.*)\n', line)
				if matchHeader:
					header = matchHeader.group(1)
					matchUniprot = re.search(r'([A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]|[O-Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9])', header)
					if matchUniprot:
						accession = matchUniprot.group(1)
						seqHeaders[accession] = header
						msaSanitized.write('>' + accession + ' HYPNO input\n')
					else:
						debugFh.write('\tMSA header at line ' + lineCount + ' missing valid Uniprot accession.\n')
						print '** HYPNO input error: given MSA headers must contain valid UniProt accessions.', \
							  '\tAccession missing for given line: ' + header
						sys.exit(1)
				else:
					msaSanitized.write(line)

	#Sanitize Tree
	if tree:
		objTree = Phylo.read(tree, "newick")

		for leaf in objTree.get_terminals():
			matchUniprot = re.search(r'([A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]|[O-Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9])', leaf.name)
			if matchUniprot:
				accession = matchUniprot.group(1)
				leafLabels[accession] = str(leaf.name)
				leaf.name = accession
			else:
				debugFh.write('\tTerminal \'' + leaf.name + '\' in the newick tree is missing a valid Uniprot accession.\n')
				print '** HYPNO input error: \'' + leaf.name + '\' in the newick tree is not a valid Uniprot accession.'
				sys.exit(1)

		Phylo.write(objTree, ID + '/' + tree, "newick")

	return leafLabels, seqHeaders

# Parses out the filename prefix 
# for the Outputs description
def parsePrefix(msaName):
	return msaName.split('.')[0]

#Calls kerf as a class to split trees into
#subtrees based on given threshold
def makeSubtrees(ID, msa, tree, threshold):
	fileMSA = str(msa)
	fileTree = str(tree)
	outputDir = str(ID)
	kerfTree = Kerf()
	kerfTree.kerfRun( fileTree, fileMSA, threshold, outputDir)

#Retreives DNA sequences for each subtree
#from Kerf generated csv MSA map file
def getDNASeqs(ID, msa, pred, execMin, seqHeaders):
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

	readProt = DNA2ProtAlign()
	ProtDict = readProt.getAAseqs(msa)

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
				nucFetchTuple = nucFetch.getSeqID(uniprotID, ProtDict[uniprotID])
				if nucFetchTuple[4] == 'null':
					missed += 1; missedList.append(uniprotID)
					if nucFetchTuple[10] != 'null':
						print '\t** HYPNO obsolete UniProt entry warning: ', nucFetchTuple[10]
						debugFh.write('\t** HYPNO obsolete UniProt entry warning: '+nucFetchTuple[10])
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
				DNASeqs2Write.append(SeqRecord(Seq(nucFetchTuple[9]), id=listLongID[i][j], description=seqHeaders[uniprotID]))
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
def alignDNASeqs(ID, msa, numSubTrees):
	for i in xrange(1, numSubTrees+1):
		fileProt = ID + '/' + msa.split('.')[0] + 'sf' + str(i) + '.' + msa.split('.')[-1]
		fileDNA = ID + '/DNAseqs' + str(i) + '.fasta'
		fileOutput = ID + '/' + 'subtree' + str(i) + '.' + msa.split('.')[-1]

		myDNA = DNA2ProtAlign()
		myDNA.alignDNAseqs(fileProt, fileDNA, fileOutput)

#For all subtrees, takes DNA alignment and
#generates new subtree using GTR algorithm
def reestimateSubtrees(ID, numSubTrees, msa):
	for i in xrange(1,numSubTrees+1):
		MSA = ID + '/' + 'subtree' + str(i) + '.' + msa.split('.')[-1]
		outputName = ID + '/' + 'subtree' + str(i) + '.ml'
		myTree = GenTree()
		myTree.makeTree(MSA , outputName, ID)


#Takes newick format output of each subtree
#and places back into site where it was
#pruned in original tree
def mergeTree(ID, tree, treeHierarchy, listLongID, MSA):
	tree = ID + '/' + tree

	myTree = GenTree()
	prunedTree = myTree.pruneTree(tree, treeHierarchy, listLongID)
	mergedTree = myTree.insertSubTrees(ID, prunedTree, treeHierarchy)

	match = re.search(r'^([^\.]+)', MSA)
	output_handle = open(match.group(1) + '.hypno.tree', "w")
	output_handle.write(mergedTree)
	output_handle.close()


def calcBranchLengths(ID, MSA, MSAtype):
	match = re.search(r'^([^\.]+)', MSA)
	tree = match.group(1) + '.hypno.tree'
	if (MSAtype == 'nuc'):
		shutil.copy(match.group(1) + '.hypno.msa', ID)
		MSA = ID + '/' + match.group(1) + '.hypno.msa'
	else:
		MSA = ID + '/' + MSA
	outputName = match.group(1) + '.opl.hypno.tree'

	myTree = GenTree()
	myTree.makeTreeBranchLengths(tree, MSA , outputName, ID, MSAtype)


def mergeAlignments(ID, MSA, numSubTrees):
	match = re.search(r'^([^\.]+)', MSA)
	with open(match.group(1) + '.hypno.msa', 'w') as outfile:
		for i in xrange(1, numSubTrees+1):
			with open(ID + '/' + 'subtree' + str(i) + '.' + MSA.split('.')[-1]) as infile:
				for line in infile:
					outfile.write(line)


def makeSingleSubtree(ID, MSA):
	shutil.copy(MSA, ID + '/' + MSA.split('.')[0] + 'sf1.' + MSA.split('.')[-1])
	with open(ID + '/' + MSA.split('.')[0] + '.csv', 'w') as outfile:
		with open(ID + '/' + MSA) as infile:
			for line in infile:
				match = re.search(r'^>.*([A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]|[O-Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9])', line)
				if match:
					outfile.write('1, ' + match.group(1) + '\n')


# Checks for proper input formatting
def validateInputs(msa, tree=None):
	# Check for existence and proper FASTA formatting of input MSA
	try:
		msaHandle = open(msa, "rU")
	except:
		print '** HYPNO input error: Given MSA file location does not exist or is not accessible: '+msa
		sys.exit(1)
	try:
		AlignIO.parse(msaHandle, "fasta").next()
	except:
		print '** HYPNO input error: improper MSA file format, must be aligned FASTA or a2m format: '+msa
		sys.exit(1)	

	if tree:
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

	if not internet_connected():
		print '** HYPNO connection error: Please connect to the internet to enable HYPNO remote database queries.'
		sys.exit(1)

	return 0

def labelLeaves(leafLabels, filename):
	treeLines = []
	with open(filename, 'r') as treeFh:
		for line in treeFh:
			treeLines.append(line)
	tree = treeLines[0]
	UniProtRegex = r'([A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]|[O-Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9])'
	for accession in re.findall(UniProtRegex, tree):
		tree = tree.replace(accession, leafLabels[accession])
	with open(filename, 'w') as treeFh:
		treeFh.write(tree)
	return 0

def internet_connected():
    try:
        response=urllib2.urlopen('http://www.google.com',timeout=1)
        return True
    except urllib2.URLError as err: pass
    return False

def main():
	#Parse arguments specifying MSA and TREE files
	parser = argparse.ArgumentParser(description='Method for HYbrid Protein NucleOtide phylogenetic gene tree reconstruction')
	parser.add_argument('--msa', type = str, help='Input Multiple Sequence Alignment (aligned FASTA or UCSC a2m format): the user must provide this argument followed by the path to an alignment of amino acid sequences, with UniProt accessions included in the sequence headers. HYPNO will retrieve the corresponding nucleic acid sequences and use the provided alignment as a template for constructing a nucleotide MSA.', required=True)
	parser.add_argument('--tree', type = str, help='Tree re-estimation (Newick format): specifying this argument along with an input tree filename results in re-estimatation of tree topology based on nucleotide sequences. If this argument is not provided, HYPNO retrieves nucleic acid sequences, overlays them on the original protein MSA and outputs this nucleotide MSA without generating a tree.', required=False)
	parser.add_argument('--k', default = 90.0, type = float, help='Subtree selection: the "subtree selection" fractional sequence identity cutoff as input to the Kerf algorithm (default: 90). To override the default setting, type "--k <X>", where X is a real number between 0 and 100. For example, to set the subtree selection value to 93, one should type "--k 93.0". (default: 90.0)', required=False)
	parser.add_argument('--n', default = 95.0, type = float, help='PercentID match: the minimum fractional sequence identity that the translation of a retrieved nucleotide sequence must have relative to the expected protein sequence for it to be accepted (default: 95). To override the default setting, type "--n <X>", where X is a real number between 0 and 100. For example, to set the percentID match value to 80, one should type "--n 80.0". (default: 95.0)', required=False)
	parser.add_argument('--s', default = 100.0, type = float, help='Retrieval minimum: the minimum fraction of nucleotide sequences that must be retrieved in order for HYPNO to continue (default: 100). In the case that a retrieval minimum of less than 100 percent is specified and HYPNO encounters a protein sequence for which a nucleotide sequence with sufficient sequence identity is unavailable, HYPNO continues executing and excludes the sequence(s) from the final nucleotide MSA and tree outputs. To override the default setting, type "--s <X>", where X is a real number between 0 and 100. For example, to set the retrieval minimum value to 90, one should type "--s 90.0". (default: 100.0)', required=False)
	parser.add_argument('--opl', default = False, action = 'store_true', help='Branch length optimization: when this argument is specified, HYPNO outputs the expected nucleotide tree as "foo.hypno.tree" but also performs midpoint rooting on the tree and calculates branch lengths using the original protein MSA, holding the tree topology fixed. The output of this branch length optimization can be found in "foo.opl.hypno.tree". (default: False)', required=False)
	parser.add_argument('--oplnuc', default = False, action = 'store_true', help='Branch length optimization: when this argument is specified, HYPNO outputs the expected nucleotide tree as "foo.hypno.tree" but also performs midpoint rooting on the tree and calculates branch lengths using the nucleotide MSA, holding the tree topology fixed. The output of this branch length optimization can be found in "foo.opl.hypno.tree". (default: False)', required=False)
	args = parser.parse_args()
	msa, tree, kerf, pred, execMin, opl, oplnuc = args.msa, args.tree, args.k, args.n, args.s, args.opl, args.oplnuc
	filePrefix = parsePrefix(msa)
	ID = time()									#Used to create directory to store files
	if tree:
		print 'Step 0 of 5: Input validation'
		validateInputs(msa,tree)
		leafLabels, seqHeaders = initialize(str(ID), msa, tree)		#Create dir and move input files
		print 'Step 1 of 5: Generating subtrees'
		makeSubtrees(str(ID), msa, tree, kerf)		#Run Kerf
		print 'Step 2 of 5: Retrieving DNA sequences'
		numSubTrees, treeHierarchy, listLongID = getDNASeqs(str(ID), msa, pred, execMin, seqHeaders)
		print 'Step 3 of 5: Mapping DNA sequences to given protein MSA'
		alignDNASeqs(str(ID), msa, numSubTrees)
		mergeAlignments(str(ID), msa, numSubTrees)
		print 'Step 4 of 5: Re-estimating subtree topologies'
		reestimateSubtrees(str(ID), numSubTrees, msa)
		print 'Step 5 of 5: Reinserting subtrees into gene tree topology'
		mergeTree(str(ID), tree, treeHierarchy, listLongID, msa)
		if opl:
			print 'Optional Step: Recalculating tree branch lengths'
			calcBranchLengths(str(ID), msa, 'prot')
		elif oplnuc:
			print 'Optional Step: Recalculating tree branch lengths'
			calcBranchLengths(str(ID), msa, 'nuc')
		labelLeaves(leafLabels, filePrefix+'.hypno.tree')
		if opl or oplnuc:
			labelLeaves(leafLabels, filePrefix+'.opl.hypno.tree')

		print '\nOutput files created in current directory:'
		if opl or oplnuc:
			print '\tHYPNO Reoptimized Hybrid Tree file: ', filePrefix+'.opl.hypno.tree'
		print '\tHYPNO Hybrid Tree file:', filePrefix+'.hypno.tree'
		print '\tNucleotide MSA file:', filePrefix+'.hypno.msa'
		print '\tHYPNO debug directory:', str(ID)
		print '\t\tDescription: please refer to HYPNO "debug" directory FAQ at'
		print '\t\thttp://tinyurl.com/bjh92an for a description of how to use'
		print '\t\tthe files provided in this directory.'
		print '\nHYPNO execution completed.'
	else:
		print 'Step 0 of 2: Input validation'
		validateInputs(msa)
		leafLabels_empty, seqHeaders = initialize(str(ID), msa)		#Create dir and move input files
		makeSingleSubtree(str(ID), msa)				#Create a CSV with entire MSA in one "subtree"
		print 'Step 1 of 2: Retrieving DNA sequences'
		numSubTrees, treeHierarchy, listLongID = getDNASeqs(str(ID), msa, pred, execMin, seqHeaders)
		print 'Step 2 of 2: Mapping DNA sequences to given protein MSA'
		alignDNASeqs(str(ID), msa, numSubTrees)
		mergeAlignments(str(ID), msa, numSubTrees)
		
		print '\nOutput files created in current directory:'
		print '\tNucleotide MSA file:', filePrefix+'.hypno.msa'
		print '\tHYPNO debug directory:', str(ID)
		print '\t\tDescription: please refer to HYPNO "debug" directory FAQ at'
		print '\t\thttp://tinyurl.com/bjh92an for a description of how to use'
		print '\t\tthe files provided in this directory.'
		print '\nHYPNO execution completed.'

if __name__ == '__main__':
	main()