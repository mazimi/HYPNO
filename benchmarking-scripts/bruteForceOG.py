#!/usr/bin/python -tt

from xml.etree.ElementTree import ElementTree
import urllib,urllib2, sys, re, xml, os
from Bio import AlignIO


#Get sequence from Uniprot given 'accession'
def getSeq(accession):
	try:
		ns = '{http://uniprot.org/uniprot}'
		pageXML = urllib.urlopen('http://www.uniprot.org/uniprot/'+accession+'.xml').read()
		tree = xml.etree.ElementTree.fromstring(pageXML)
		items = tree.getiterator(ns+'uniprot')
		item = items[0]
		sequence = item.find(ns+'entry').find(ns+'sequence')
	except:
		return 'ErrNoSeq'

	return sequence.text.strip()


#Get PWID of two sequences (copied from RTfetch.py)
def getpID(mySeq1,mySeq2):											
			seqA = open('tempA.fasta', 'w')			# Print sequences to temporary files for EMBOSS needle call
			seqB = open('tempB.fasta', 'w')				
			print >>seqA, '> seqA\n'+mySeq1+'\n'
			print >>seqB, '> seqB\n'+mySeq2+'\n'
			seqA.close()							# Close temporary file handle references
			seqB.close()
			os.system('needle -asequence tempA.fasta -sprotein1 -bsequence tempB.fasta -sprotein2 -gapopen 10 -gapextend 0.5 -outfile tempAlign.needle -auto')
			needle = open('tempAlign.needle','rU')
			alignment = AlignIO.read(needle,"emboss")# AlignIO BioPython Module reads out EMBOSS globally aligned sequences
			i=0 									 # Global alignment --> only 1 counter necessary for both sequences
			counter = 0 							 # Global counter
			gaps = 0
			sequence0 = alignment[0]
			sequence1 = alignment[1]
			seq0 = str(sequence0.seq)
			seq1 = str(sequence1.seq)
			list0 = list(seq0)
			list1 = list(seq1)
			while counter < len(list0):
				topAA = list0[counter]
				bottomAA = list1[counter]
				# Considers gaps and mismatches in both sequences for computing percent identity
				if topAA != bottomAA:
					gaps = gaps + 1
					pass
				else:
					i = i + 1
				counter = counter + 1
			percent = 100*i/len(seq0)
			os.remove('tempA.fasta')	# Remove these temporary files that you no longer need
			os.remove('tempB.fasta')
			return percent


#Determine indices in matrixPWID that correspond to accession combination
def buildGroupScore(groups):
	listAccCombos = []
	scoreIndices = []

	if len(groups) > 2:
		(subs, subScores) = buildGroupScore(groups[1:])
	else:
		(subs, subScores) = [groups[1], range(len(groups[1]))]

	for acc in groups[0]:
		for sub in subs:
			listAccCombos.append(acc + ',' + sub)

	for i,acc in enumerate(groups[0]):
		for subInd in subScores:
			scoreIndices.append(str(i) + ',' + str(subInd))

	return (listAccCombos, scoreIndices)


# Return OG with highest total PWID between all accessions
def getBestOG(strAllAcc):
	numTaxa = 0
	matrixAccessions = []
	matrixSeqs = []
	matrixPWID = []
	listMaxAccessions = []

	# Build matrix of accesions
	matrixAccessions = eval(strAllAcc)
	#DEBUG: Print accession matrix
	#print matrixAccessions

	# Build matrix of sequences
	for taxa in matrixAccessions:
		listSeqs = []
		for accession in taxa:
			mySeq = getSeq(accession)
			if mySeq == 'ErrNoSeq':
				print 'Error: Couldn\'t find sequence for accession: ' + accession
				listSeqs.append("XXX")
			else:
				listSeqs.append(mySeq)
		matrixSeqs.append(listSeqs)

	# Build matrix of all percent pairwise identity
	matrixPWID = [[[] for i in range(len(matrixSeqs))] for j in range(len(matrixSeqs))]
	for i in xrange(0,len(matrixSeqs)):
		for j in xrange(i+1,len(matrixSeqs)):
			matrixSeqs1Seq2 = []
			for seq1 in matrixSeqs[i]:
				listSeqs2 = []
				for seq2 in matrixSeqs[j]:
					if (seq1 == "XXX") or (seq2 == "XXX"):
						listSeqs2.append(0)
					else:
						listSeqs2.append(getpID(seq1,seq2))
				matrixSeqs1Seq2.append(listSeqs2)
			matrixPWID[i][j] = matrixSeqs1Seq2


	#DEBUG Print PWID matrix
	#for i in xrange(0,len(matrixSeqs)):
	#	for j in xrange(0,len(matrixSeqs)):
	#		print str(matrixPWID[i][j]) + '\t\t',
	#	print '\n'

	# Get all combinations of proteins and indices of PWID for combination
	(listAccCombos,scoreIndices) = buildGroupScore(matrixAccessions)

	# Compute sum of PWID for each accession combination and select highest scoring
	maxScore = 0
	maxIndex = 0
	for m,indices in enumerate(scoreIndices):
		listIndices = indices.split(',')
		currSum = 0
		minScore = 100
		for i in xrange(0,len(matrixSeqs)):
			for j in xrange(i+1,len(matrixSeqs)):
				currScore = matrixPWID[i][j][eval(listIndices[i])][eval(listIndices[j])]
				currSum += currScore
				if currScore < minScore:
					minScore = currScore

		#DEBUG: Print all combinations and sum of PWIDs
		#print listAccCombos[m], currSum, minScore
		if currSum > maxScore:
			maxScore = currSum
			maxIndex = m
			minPID = minScore


	return (listAccCombos[m].split(','), minPID)


def main(inFile,outFile):
	allAccessions = []

	with open(inFile, 'rU') as fin:
		for line in fin:
			match = re.search(r'^(\S+)\t(\S+)\t(\S+)\t(.+)\n$', line)
	
			print match.group(4)
			(bestOG,minPID) = getBestOG(match.group(4))
	
			flagDuplicate = False
			for accession1 in bestOG:
				for accession2 in allAccessions:
					if accession1 == accession2:
						flagDuplicate = True
	
			if flagDuplicate == False:
				with open(outFile, 'a') as foutAccept:
					foutAccept.write(match.group(1) + '\t' + match.group(2) + '\t' + match.group(3) + '\t' + str(bestOG) + '\t' + str(minPID) + '\n')
					allAccessions = allAccessions + bestOG
			else:
				with open(outFile + '.redundant', 'a') as foutReject:
					foutReject.write(match.group(1) + '\t' + match.group(2) + '\t' + match.group(3) + '\t' + str(bestOG) + '\t' + str(minPID) + '\n')

if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2])
