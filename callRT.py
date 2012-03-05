#!/usr/bin/python -tt
# Copyright Nima Emami, 2012

from RTfetch import RTfetch
import sys

# main: Delegates nucleotide sequence fetch
# input: none
# output: none
# NOTE: Programmatic Uniprot Access: http://www.uniprot.org/faq/28#id_mapping_examples
def main():
	uniprotIDs = ['Q6GNU6', 'E1BR73', 'E1BS49', 'D2HE46', 'A7E2A2']
	nucFetch = RTfetch()
	for uniprotID in uniprotIDs:
		print ''
		dataTuple = nucFetch.getSeqID(uniprotID)
		print dataTuple

if __name__ == '__main__':
  main()