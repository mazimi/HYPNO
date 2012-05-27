#!/opt/python/bin/python2.7

#Script that takes in a protein Multiple Sequence Alignment - MSA - (FASTA format)
#and a file containing the corresponding DNA sequences (FASTA format) and reverse-
#translates the protein to DNA, preserving indel characters and overall alignment.

from Bio import SeqIO, Seq, Phylo
import sys

class DNA2ProtAlign:

	# __init__: Constructor method
	# input: self
	# output: none
	def __init__(self):
		pass

	def alignDNAseqs(self, fileProt, fileDNA, fileOutput, fileAllTree):

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

		#Retreive entire DNA record from file to use for FASTA reconstruction
		#record.seq will be overwritten, other information will be preserved
		def getDNArecord(fileDNA):
			seqDNAid = []
			handle = open(fileDNA, "rU")
			for record in SeqIO.parse(handle, "fasta"):
				seqDNAid.append(record)
			handle.close()

			return seqDNAid

		#Reverse-translate protein to corresponding DNA sequence from retreived
		#DNA sequence.
		def translateProt2DNAmsa(seqProt, seqDNA, idProt):
			checkCodons = False	#Set 'True' to check that codons match residues
			gapHandlingDNA = 2	#Set '1' to transfer gap to protein
								#Set '2' to use aligned codon from nearest neighbor
								#Set '3' TBD

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
						#TODO: How should mismatches be handled?
						if str(seqDNA[i][indexDNA*3:indexDNA*3+3]) == '***':
							if gapHandlingDNA == 1:		#Transfer DNA gap to protein
								currSeq += '---'
							if gapHandlingDNA == 2:
								fullTree = Phylo.read(fileAllTree, 'newick')
								#DEBUG: Draw an ASCII tree, just for Nima
								print Phylo.draw_ascii(fullTree)
								minDistance = 99999					#Set very high value
								nearestNeighbor = -1				#Start at -1 in case no neighbor found
								#Loop over all proteins in clade looking for nearest neighbor with matching AA at the
								#same position and transfer codon if available.
								for k in xrange(0, len(seqProt)):
									if (idProt[i] != idProt[k]) and (seqProt[k][j] == seqProt[i][j] ):	#TODO: Currently only checking for AA w/ same case
										distNodes = fullTree.distance(idProt[i], idProt[k])
										#DEBUG: Print distance values
										print "Distance between " + str(idProt[i]) + " and " + str(idProt[k])
										print distNodes
										if distNodes < minDistance:		#Select current leaf it's the closer to target than last nearest
											minDistance = distNodes
											nearestNeighbor = k
								if nearestNeighbor > -1:		#A nearest neighbor with matching AA was found!
									#DEBUG: List nearest neighbor with matching AA
									print "Nearest neighbor is: " + str(idProt[nearestNeighbor])
									currSeq += seqDNA[nearestNeighbor][indexDNA*3:indexDNA*3+3]
						else:
							if checkCodons:
								if str(Seq.translate(seqDNA[i][indexDNA*3:indexDNA*3+3])) == str(seqProt[i][j]):
									currSeq += seqDNA[i][indexDNA*3:indexDNA*3+3]
								else:
									currSeq += '!!!'
							else:
								currSeq += seqDNA[i][indexDNA*3:indexDNA*3+3]
						indexDNA += 1
					else:
						#Grab next Lower Case codon in seq, check to make sure it codes for current AA (if checkCodons=True),
						#if match, place Lower Case codon in DNA MSA sequence, if no match, place '!!!',
						#move on to next codon.
						#TODO: How should mismatches be handled?
						if str(seqDNA[i][indexDNA*3:indexDNA*3+3]) == '***':
							if gapHandlingDNA == 1:		#Transfer DNA gap to protein
								currSeq += '---'
							if gapHandlingDNA == 2:
								fullTree = Phylo.read('bpg0233328.ml', 'newick')
						else:
							if checkCodons:
								upperSeqProtein = seqProt[i][j].upper()
								if str(Seq.translate(seqDNA[i][indexDNA*3:indexDNA*3+3])) == str(upperSeqProtein):
									currSeq += seqDNA[i][indexDNA*3:indexDNA*3+3].lower()
								else:
									currSeq += '!!!'
							else:
								currSeq += seqDNA[i][indexDNA*3:indexDNA*3+3].lower()
						indexDNA += 1
				msaDNA.append(currSeq)

			return msaDNA


		seqProt = getProtSeq(fileProt)
		idProt = getProtID(fileProt)
		seqDNA = getDNASeq(fileDNA)
		seqDNArecord = getDNArecord(fileDNA)

		msaDNA = translateProt2DNAmsa(seqProt, seqDNA, idProt)

		if len(msaDNA) != len(seqProt):
			print "Number of protein and DNA sequences does not match!"
		else:
			#TODO: output to file/return
			for i in xrange(0, len(msaDNA)):
				seqDNArecord[i].seq = msaDNA[i]
			output_handle = open(fileOutput, "w")
			SeqIO.write(seqDNArecord, output_handle, "fasta")
			output_handle.close()

if __name__ == '__main__':
	sys.exit(1)