#!/opt/python/bin/python2.7

#AUTHORS: mazimi, nemami. 2012-6-10
#LICENSE: https://github.com/mazimi/HYPNO/blob/master/LICENSE
#DESCRIPTION: Parser for OrthoDB flatfile to extract Orthology Group (OG)
# for each UniProt Accession (UniProt_ACC) in tab delimited format.
#LATEST VERSION: https://github.com/mazimi/HYPNO/blob/master/OrthoGroups/OrthoParse.py

import argparse
import re

def main():
	#Parse arguments specifying input and output files
	parser = argparse.ArgumentParser(description='Parser for OrthoDB flatfile to extract Orthology Group (OG) for each UniProt Accession (UniProt_ACC) in tab delimited format.')
	parser.add_argument('--input', required=True)
	parser.add_argument('--output', required=True)
	args = parser.parse_args()
	inFile = str(args.input)
	outFile = str(args.output)

	#Open files for input/output/error
	fin = open(inFile, 'rU')
	fout = open(outFile, 'w')
	ferr = open('OrthoParse.error', 'w')

	#Read input file line by line. Check that each line contains the correct number of entries,
	# capturing Orthology Group ID and UniProt Accession ID. UniProt Accession IDs are checked
	# for validity. If everything checks out, entry is written to output file. Incorrect number
	# of entries in a line or invalid UniProt Accession IDs are tracked in the error output file.
	numLine = 1
	for line in fin:
		#OrthoDB flatfiles contain 7 tab delimited entries (see first line of file for labels).
		match = re.search(r'[^\t]*\t([^\t]*)\t[^\t]*\t[^\t]*\t[^\t]*\t([^\t]*)\t[^\t]*\n', line)
		if match:
			#Check for valid UniProt Accession ID, see here for guidelines: http://www.uniprot.org/manual/accession_numbers
			matchUniprot = re.search(r'^([A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]|[O-Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9])', match.group(2))
			if matchUniprot:
				fout.write(match.group(2) + '\t' + match.group(1) + '\n')
			else:
				ferr.write('UniProt ID Error: \'' + match.group(2) + '\' on line: ' + str(numLine) + '\n')
		else:
			ferr.write('Parse Error: line ' + str(numLine) + '\n')
		numLine += 1

	#Close files after input/output/error
	fin.close()
	fout.close()
	ferr.close()

if __name__ == '__main__':
	main()