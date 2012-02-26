#!/opt/python/bin/python2.7

#Script that takes in a protein Multiple Sequence Alignment - MSA - (FASTA format)
#and a file containing the corresponding DNA sequences (FASTA format) and reverse-
#translates the protein to DNA, preserving indel characters and overall alignment.

import argparse
from Bio import SeqIO, Seq
import sys


#Read protein sequences from file
def getProtSeq(fileProt):
	seqProt = []
	handle = open(fileProt, "rU")
	for record in SeqIO.parse(handle, "fasta"):
		seqProt.append(record.seq)
	handle.close()

	return seqProt

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

#Reverse-translate protein to corresponding DNA sequence using the codon
#that is found to occur in nature.
def translateProt2DNAmsa(seqProt, seqDNA):
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
			else:
				#Grab next codon in seq, check to make sure it codes for current AA,
				#if match, place codon in DNA MSA sequence, if no match, place '!!!',
				#move on to next codon.
				#TODO: How should mismatches be handled?
				if str(Seq.translate(seqDNA[i][indexDNA*3:indexDNA*3+3])) == str(seqProt[i][j]):
					currSeq += seqDNA[i][indexDNA*3:indexDNA*3+3]
				else:
					currSeq += '!!!'
				indexDNA += 1
		msaDNA.append(currSeq)

	return msaDNA

def main():	
	#Parse arguments specifying protein MSA file (FASTA) and DNA sequence file (FASTA)
	#TODO: This script will become a standalone module.
	parser = argparse.ArgumentParser(description='Method for translating protein MSA to DNA MSA')
	parser.add_argument('--prot-file', required=True)
	parser.add_argument('--DNA-file', required=True)
	parser.add_argument('--output-file', required=True)
	args = parser.parse_args()
	fileProt = str(args.prot_file)
	fileDNA = str(args.DNA_file)
	fileOutput = str(args.output_file)

	seqProt = getProtSeq(fileProt)
	seqDNA = getDNASeq(fileDNA)
	seqDNArecord = getDNArecord(fileDNA)

	msaDNA = translateProt2DNAmsa(seqProt, seqDNA)

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
	main()