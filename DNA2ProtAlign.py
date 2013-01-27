#!/usr/bin/env python2.7

#Script that takes in a protein Multiple Sequence Alignment - MSA - (FASTA format)
#and a file containing the corresponding DNA sequences (FASTA format) and reverse-
#translates the protein to DNA, preserving indel characters and overall alignment.

from Bio import SeqIO, Seq, Phylo
import sys, re

class DNA2ProtAlign:

	# __init__: Constructor method
	# input: self
	# output: none
	def __init__(self):
		pass

	def getAAseqs(self, fileProt):

				#Read protein sequences from file
		def getProtDict(fileProt):
			ProtDict = {}
			with open(fileProt, "rU") as handle:
				for record in SeqIO.parse(handle, "fasta"):
					accession = re.search(r'([A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]|[O-Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9])', record.id)
					if accession:
						ProtDict[accession.group(1)] = str(record.seq).replace('-','')

			return ProtDict

		return getProtDict(fileProt)

	def alignDNAseqs(self, fileProt, fileDNA, fileOutput):

		#Read protein sequences from file
		def getProtSeq(fileProt):
			seqProt = []
			handle = open(fileProt, "rU")
			for record in SeqIO.parse(handle, "fasta"):
				seqProt.append(record.seq)
			handle.close()

			return seqProt

		#Read protein IDs from file
		def getProtID(fileProt):
			seqID = []
			handle = open(fileProt, "rU")
			for record in SeqIO.parse(handle, "fasta"):
				seqID.append(record.id)
			handle.close()

			return seqID

		#Read DNA sequences from file
		def getDNASeq(fileDNA):
			seqDNA = []
			handle = open(fileDNA, "rU")
			for record in SeqIO.parse(handle, "fasta"):
				seqDNA.append(record.seq)
			handle.close()

			return seqDNA

		#Read DNA IDs from file
		def getDNAid(fileDNA):
			seqDNAid = []
			handle = open(fileDNA, "rU")
			for record in SeqIO.parse(handle, "fasta"):
				seqDNAid.append(record.id)
			handle.close()

			return seqDNAid

		#Retreive entire DNA record from file to use for FASTA reconstruction
		#record.seq will be overwritten, other information will be preserved
		def getDNArecord(fileDNA):
			seqDNAid = []
			handle = open(fileDNA, "rU")
			for record in SeqIO.parse(handle, "fasta"):
				seqDNAid.append(record)
			handle.close()

			return seqDNAid

		def filterProtSeqs(idProt, idDNA, seqProt):
			seqProtFiltered = []
			idProtFiltered = []
			for i in xrange(0,len(idDNA)):
				matchDNA = re.search(r'([A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]|[O-Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9])', idDNA[i])
				if matchDNA:
					for j in xrange(0,len(idProt)):
						matchProt = re.search(r'([A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]|[O-Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9])', idProt[j])
						if matchProt:
							if matchDNA.group(1) == matchProt.group(1):
								seqProtFiltered.append(seqProt[j])
								idProtFiltered.append(matchProt.group(1))
								break

			return seqProtFiltered, idProtFiltered

		#Reverse-translate protein to corresponding DNA sequence from retreived
		#DNA sequence.
		def translateProt2DNAmsa(seqProt, seqDNA, idProt):
			checkCodons = False	#Set 'True' to check that codons match residues

			#Make sure number of DNA sequences matches number of proteins aligned
			if len(seqProt) != len(seqDNA):
				return

			msaDNA = []

			for i in xrange(0, len(seqProt)):
				currSeq = ''
				indexDNA = 0
				for j in xrange(0, len(seqProt[i])):
					if seqProt[i][j] == '-':
						currSeq += '---'
					elif seqProt[i][j] =='.':
						currSeq += '...'
					elif seqProt[i][j].isupper() == True:
						#Grab next Upper Case codon in seq, check to make sure it codes for current AA (if checkCodons=True),
						#if match, place Upper Casecodon in DNA MSA sequence, if no match, place '!!!',
						#move on to next codon.
						if str(seqDNA[i][indexDNA*3:indexDNA*3+3]) == '***':
							currSeq += '---'
						else:
							currSeq += seqDNA[i][indexDNA*3:indexDNA*3+3]
						indexDNA += 1
					else:
						#Grab next Lower Case codon in seq, check to make sure it codes for current AA (if checkCodons=True),
						#if match, place Lower Case codon in DNA MSA sequence, if no match, place '!!!',
						#move on to next codon.
						if str(seqDNA[i][indexDNA*3:indexDNA*3+3]) == '***':
							currSeq += '---'
						else:
							currSeq += seqDNA[i][indexDNA*3:indexDNA*3+3].lower()
						indexDNA += 1
				msaDNA.append(currSeq)

			return msaDNA


		seqProt = getProtSeq(fileProt)
		idProt = getProtID(fileProt)
		seqDNA = getDNASeq(fileDNA)
		idDNA = getDNAid(fileDNA)
		seqDNArecord = getDNArecord(fileDNA)
		seqProtFiltered, idProtFiltered = filterProtSeqs(idProt, idDNA, seqProt)


		msaDNA = translateProt2DNAmsa(seqProtFiltered, seqDNA, idProtFiltered)

		if len(msaDNA) != len(seqProtFiltered):
			print "Number of protein and DNA sequences does not match!"
		else:
			for i in xrange(0, len(msaDNA)):
				seqDNArecord[i].seq = msaDNA[i]
			output_handle = open(fileOutput, "w")
			SeqIO.write(seqDNArecord, output_handle, "fasta")
			output_handle.close()

if __name__ == '__main__':
	sys.exit(1)