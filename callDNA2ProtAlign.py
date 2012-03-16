#!/opt/python/bin/python2.7

from DNA2ProtAlign import DNA2ProtAlign
import sys

# main: Assign protein MSA (a2m), DNA sequences (FASTA) and aligned DNA output file (FASTA)
# input: none
# output: none
def main():
	fileProt = "protAlign.a2m"
	fileDNA = "DNAseqs.fasta"
	fileOutput = "DNAaligned.fasta"
	myDNA = DNA2ProtAlign()
	myDNA.alignDNAseqs(fileProt, fileDNA, fileOutput)

if __name__ == '__main__':
  main()