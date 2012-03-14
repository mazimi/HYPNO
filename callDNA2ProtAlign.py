#!/opt/python/bin/python2.7

from DNA2ProtAlign import DNA2ProtAlign
import sys

# main: Delegates nucleotide sequence fetch
# input: none
# output: none
# NOTE: Programmatic Uniprot Access: http://www.uniprot.org/faq/28#id_mapping_examples
def main():
	fileProt = "protAlign.a2m"
	fileDNA = "DNAseqs.fasta"
	fileOutput = "DNAaligned.fasta"
	myDNA = DNA2ProtAlign()
	myDNA.alignDNAseqs(fileProt, fileDNA, fileOutput)

if __name__ == '__main__':
  main()