#!/opt/python/bin/python2.7

from DNA2ProtAlign import DNA2ProtAlign
import sys

# main: Assign protein MSA (a2m), DNA sequences (FASTA) and aligned DNA output file (FASTA)
# input: none
# output: none
def main():
	fileProt = "bpg0233328sf3.a2m"
	fileDNA = "DNAseqs3.fasta"
	fileOutput = "DNAaligned3.fasta"
	fileAllMSA = "bpg0233328.a2m"
	fileAllTree = "bpg0233328.ml"
	myDNA = DNA2ProtAlign()
	myDNA.alignDNAseqs(fileProt, fileDNA, fileOutput, fileAllTree)

if __name__ == '__main__':
  main()