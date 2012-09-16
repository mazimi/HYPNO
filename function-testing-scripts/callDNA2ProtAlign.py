#!/opt/python/bin/python2.7

from DNA2ProtAlign import DNA2ProtAlign
import sys

# main: Assign protein MSA (a2m), DNA sequences (FASTA) and aligned DNA output file (FASTA)
# input: none
# output: none
def main():
	fileProt = "1333735262.28 copy/bpg0240116sf4.a2m"
	fileDNA = "1333735262.28 copy/DNAseqs4.fasta"
	fileOutput = "1333735262.28 copy/DNAaligned4.fasta"
	myDNA = DNA2ProtAlign()
	myDNA.alignDNAseqs(fileProt, fileDNA, fileOutput)

if __name__ == '__main__':
  main()