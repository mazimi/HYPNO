#!/usr/bin/python -tt
# Copyright Nima Emami, 2012

from RTfetch import RTfetch
import sys, time

# main: Delegates nucleotide sequence fetch
# input: none
# output: none
# NOTE: Programmatic Uniprot Access: http://www.uniprot.org/faq/28#id_mapping_examples
def main(infile, outfile):
	lineList = []
	fhIn = open(infile, 'rU')
	for line in fhIn:
		lineList.append(line)
	fhIn.close()
	firstline = lineList.pop(0)
	uniprotIDs = firstline.split()
	nucFetch = RTfetch()
	fh = open(outfile, 'w')
	start_time = time.time()
	for uniprotID in uniprotIDs:
		start = time.time()
		print >>fh,""
		print uniprotID
		dataTuple = nucFetch.getSeqID(uniprotID)
		print >>fh, dataTuple
		end = time.time()-start
		print >>fh, "elapsed: "+str(end)
	end_time = time.time()-start_time
	print 'Total Runtime: '+str(end_time)
	fh.close()

if __name__ == '__main__':
  sys.argv.pop(0)
  infile = sys.argv.pop(0)
  outfile = sys.argv.pop(0)
  main(infile, outfile)
  