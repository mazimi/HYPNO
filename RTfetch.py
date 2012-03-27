#!/usr/bin/python -tt
# Copyright Nima Emami, 2012

from Bio import Entrez, SeqIO, AlignIO
from xml.etree.ElementTree import ElementTree
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
import urllib,urllib2, sys, re, xml, os

# NOTES: Requires ElementTree module, EMBOSS package, and BioPython for functionality

class RTfetch:

	# __init__: Constructor method, delegates nucleotide sequence fetch
	# input: self
	# output: none
	def __init__(self):
		Entrez.email = 'timanoed154@yahoo.com'				# mandatory NCBI email for reference

	# getSeqID: Delegates nucleotide sequence fetch
		# input: uniprotID (0)
		# output: prints uniprot ID, nucleotide ID + DNA ORF matches to standard output, or error message if unsuccesful
		# NOTE: Programmatic Uniprot Access: http://www.uniprot.org/faq/28#id_mapping_examples
			# Tuple indices:
			# (1) Uniprot ID query
			# (2) Nucleotide ID (default = 'null')
			# (3) Database Source (EMBL, NCBI, TBLASTN, default = 'null')
			# (4) Protein sequence
			# (5) DNA Sequence (default = 'null')
			# (6) Aligned protein hit (default = 'null')
			# (7) Percent Identity (default = 'null')
			# (8) ORF Starting Index (default = 'null')
			# (9) ORF Ending Index (default = 'null')
			# (10) mapped DNA sequence (default = 'null')
	def getSeqID(self, uniprotID):
		
		# orfFinder: 6 frame ORF search on DNA/mRNA sequence (Adapted & Modified from BioPython)
		# input: Nucleic Acid Sequence (0), Translation Table Number (1), ORF lenght lower bound (2),
		#			Protein Sequence for alignment confirmation (3)
		# output: Start Index of Hit (0), End Index of Hit (1), flag = 1 (2) OR flag = 0 (0)
		# NOTE: Source code: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc224
			# NOTE: Translation tables: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
		# May need to be modified based on taxonomic classifications, currently set to classical genetic code
		def orfFinder(seq, trans_table, min_protein_length, proSeq, source):
			bestSoFar = 0 				# Set defaults
			startSoFar = 0 				# ^
			pStartSoFar = 0 			# |
			pEndSoFar = 0 				# |
			endSoFar = 0 				# |
			pidSoFar = 0 				# |
			alignedSoFar = '' 			# |
			seq_len = len(seq)			# Length of original protein sequence, for the purpose of indexing start and end indices
			for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:	# Forward vs. Reverse frames
				for frame in range(3):											# 3 reading frames for each
					trans = str(nuc[frame:].translate(trans_table))			# ** Translation based on trans_table not working...
					trans_len = len(trans)									# --> Quick fix: check for initial methionine
					aa_start = 0
					aa_end = 0
					while aa_start < trans_len:
						aa_end = trans.find("*", aa_start)					# Find Stop codon
						if aa_end == -1:
							aa_end = trans_len
						# ** Current Heuristic: Only look for DNA sequences which code for proteins 3/4ths to 5/4ths of the size
						# -->  May want to refine this at some point
						if trans[aa_start] == 'M' and 1.25*min_protein_length >= aa_end-aa_start >= 0.75*min_protein_length:
							(percentID, aligned) = getpID(trans[aa_start:aa_end].lower(),proSeq.lower())
							if int(percentID) > 99:					# If it is almost a perfect match
								if strand == 1:						# Forward frame
									start = max(0,frame+aa_start*3)			# start index (w/ max because aa_start + frame may be less than 0)
									end = min(seq_len,frame+aa_end*3+3)		# end index (w/ +3 at end because frame goes negative instead of positive)
									mappedDNA = mapDNA(start,end,seq,aligned)	# Call DNA mapping subroutine to map DNA to aligned protein to query
								else:								# Reverse frame
									start = max(0,seq_len-frame-aa_end*3-3)	# Backwards case, so everything is effectively reversed
									end = seq_len-frame-aa_start*3 			# end index is now at the beginning, so subtract off AA start
									mappedDNA = mapDNA(start,end,seq,aligned)	# Call DNA mapping subroutine to map DNA to aligned protein to query
									mappedSeq = Seq(mappedDNA)					# RevComp the DNA sequence
									mappedSeq = mappedSeq.reverse_complement()
								return (end,start,aligned,int(percentID),str(mappedDNA)) # Flip Start and End indices because of start and stop index
							elif int(percentID) > bestSoFar:	# Trying to maximize percent identity
																# (See isoform discussion: http://www.uniprot.org/faq/30)
								bestSoFar = int(percentID)
								if strand == 1:
									start = max(0,frame+aa_start*3)
									end = min(seq_len,frame+aa_end*3+3)
								else:
									start = max(0,seq_len-frame-aa_end*3-3)
									end = seq_len-frame-aa_start*3                     
								startSoFar = start
								endSoFar = end
								alignedSoFar = aligned
								pidSoFar = percentID
						aa_start = aa_start+1
			mappedDNA = mapDNA(startSoFar,endSoFar,seq,alignedSoFar)
			return (startSoFar,endSoFar,alignedSoFar,int(pidSoFar),mappedDNA)

		# mapDNA: mapping DNA sequence to the Needleman-Wunsch aligned protein to the original protein query
		# input: start index (0), end index (1), dna (2), aligned protein (3)
		# output: mapped DNA sequence
		def mapDNA(startIndex,endIndex,dna,alignedProtein):
			dna = str(dna)
			mapped = ''
			dnaCounter = startIndex
			protCounter = 0
			while(protCounter < len(alignedProtein)):	# While you have not surpassed the length of the aligned protein sequence
				if(alignedProtein[protCounter] == ' ' or alignedProtein[protCounter] == '-'):
					mapped = mapped + '   '				# For each gap in the protein alignment, insert three gaps in the DNA alignment
				else:
					mapped = mapped + dna[dnaCounter:dnaCounter+3]	# Otherwise, insert the next three nucleic acids
					dnaCounter = dnaCounter + 3 					# and increment both counters
				protCounter = protCounter + 1 							# instead of just one.
			return mapped

		# pID: pairwise alignment subroutine (Adapted and modified from below source)
		# input: 2 nucleotide sequences (0,1)
		# output: integer value between 0 and 100 (0)
		# NOTE: Uses EMBOSS package needle executable call
		# NOTE: Prints to temporary file handle tempAlign.needle
		def getpID(target,template):											
			seqA = open('tempA.fasta', 'w')			# Print sequences to temporary files for EMBOSS needle call
			seqB = open('tempB.fasta', 'w')				
			print >>seqA, '> seqA\n'+target+'\n'
			print >>seqB, '> seqB\n'+template+'\n'
			seqA.close()							# Close temporary file handle references
			seqB.close()
			os.system('needle -asequence tempA.fasta -sprotein1 -bsequence tempB.fasta -sprotein2 -gapopen 10 -gapextend 0.5 -outfile tempAlign.needle > nul')
			needle = open('tempAlign.needle','rU')
			alignment = AlignIO.read(needle,"emboss")		# AlignIO BioPerl Module reads out EMBOSS globally aligned sequences
			i=0 									 # Global alignment --> only 1 counter necessary for both sequences
			counter = 0 							 # Global counter
			sequence0 = alignment[0]
			sequence1 = alignment[1]
			seq0 = str(sequence0.seq)
			seq1 = str(sequence1.seq)
			list0 = list(seq0)
			list1 = list(seq1)
			while counter < len(list0):
				topAA = list0[counter]
				bottomAA = list1[counter]
				# Considers gaps and mismatches in both sequences for computing percent identity
				if topAA == '-' or topAA == ' ' or  bottomAA == '-' or bottomAA == ' ' or topAA != bottomAA:
					pass
				else:
					i = i + 1
				counter = counter + 1
			percent = 100*i/len(seq0)
			os.remove('tempA.fasta')	# Remove these temporary files that you no longer need
			os.remove('tempB.fasta')
			return (str(percent),seq0.upper())

		# orfLength: determines length of ORFs based on protein sequence length via uniprot query
		# input: uniprot ID (0)
		# output: tuple of (length of protein, protein sequence)
		def orfLength(protID):
			page = urllib.urlopen('http://www.uniprot.org/uniprot/'+protID+'.fasta').read()
			pageLines = page.split('\n')
			pageLines.pop(0)						# Remove fasta description > ...
			merged = ''.join(pageLines)				
			merged = merged.strip()					# Remove whitespace
			proTup = (len(merged),merged)
			return proTup
			
		# getSeq: use NCBI accession retrieved by tBLASTn to do a nucleotide query 
		# input: NCBI accession (0)
		# output: returns tuple of (NCBI_ID,'NCBI')
		def getSeq(NCBI_ID):
			try:
				handle = Entrez.efetch(db="nucleotide", id=NCBI_ID,rettype="fasta")	# use e-fetch to get sequence
				record = SeqIO.read(handle, "fasta")		# use SeqIO to get out the fasta DNA/mRNA sequence
				Seq = record.seq
				dna_seq = Seq.tostring()
			except:
				# print >>fh,'problem reading: '+url+' or performing Entrez efetch' # catch URL connection errors
				# sys.stderr.write('problem reading: '+url+' or performing Entrez efetch')	# catch URL connection errors
				dna_seq = ''
			return dna_seq

		# taxBLASTn: tries does uniprot --> NBCI nucleotide ID fetch
		# input: uniprot ID (0), common taxonomic organism name (1)
		# output: tuple (uniprot ID, NCBI nucleotid ID,DNA ORF matches to standard output
		# NOTE: Biopython BLAST documentation: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc81
		def taxBLASTn(uniprotID,taxName,proSeq):
			mapTuple = ()
			try:
				# print >>fh,'\nFor '+uniprotID+' -- '+taxName+', performing taxBLASTn:'
				result_handle = NCBIWWW.qblast("tblastn", "nr", proSeq, expect = .0001, entrez_query = taxName+'[organism]')
				string = result_handle.read()
				result_handle.close()
				tree = xml.etree.ElementTree.fromstring(string)
				iteration = tree.find("BlastOutput_iterations/Iteration")
				hits = iteration.findall("Iteration_hits/Hit")
				topHit = hits[0]
				accessionNCBI = topHit.findtext("Hit_accession")
				qseq = topHit.findtext("Hit_hsps/Hsp/Hsp_qseq")
				hseq = topHit.findtext("Hit_hsps/Hsp/Hsp_hseq")
				midseq = topHit.findtext("Hit_hsps/Hsp/Hsp_midline")
				Hit_id = topHit.findtext("Hit_id")
				Hit_from = int(topHit.findtext("Hit_hsps/Hsp/Hsp_hit-from"))
				Hit_to = int(topHit.findtext("Hit_hsps/Hsp/Hsp_hit-to"))
				match = re.search(r'gi\|(\w+)\|',Hit_id)
				GI = match.group(1)
				# print >>fh,'Q: '+qseq
				# print >>fh,"M: "+midseq
				# print >>fh,"H: "+hseq
				sec = accessionNCBI[1:2]
				onetwo = accessionNCBI[0:2]
				# TODO: improve this, make more robust
				flagOne = sec == 'C' or sec == 'G' or sec == 'T' or sec == 'W' or sec == 'Z' or sec == 'S'
				flagTwo = onetwo == 'AP' or onetwo == 'BS' or onetwo == 'AL' or onetwo == 'BX' or onetwo == 'CR'\
							or onetwo == 'CT' or onetwo == 'CU' or onetwo == 'FP' or onetwo == 'FQ' or onetwo == 'FR'\
							or onetwo == 'AE' or onetwo == 'CP' or onetwo == 'CY' 
				if flagOne or flagTwo:
					dna_seq = chromParse(GI,Hit_from,Hit_to)
				else:
					dna_seq = getSeq(accessionNCBI)
				if dna_seq != '':
					mapTuple = (uniprotID,accessionNCBI,dna_seq,'TBLASTN',midseq)
					return mapTuple
				else:
					# print >>fh,('\tFor uniprot ID '+uniprotID+' and NCBI accession '+accessionNCBI+', problem with NCBI efetch sequence retrieval.')	# catch URL connection errors
					# sys.stderr.write('\tFor uniprot ID '+uniprotID+' and NCBI accession '+accessionNCBI+', problem with NCBI efetch sequence retrieval.')	# catch URL connection errors
					return mapTuple
			except:
				# print >>fh,'\tProblem with a) taxon access from Uniprot,\
				 				 # \n\tb) reading organism specific tBLASTn,\
								 # \n\tc) XML parsing of NCBI accession ID for the top hit' # catch URL connection errors
				# sys.stderr.write('\tProblem with a) taxon access from Uniprot,\
				# 				 \n\tb) reading organism specific tBLASTn,\
				# 				 \n\tc) XML parsing of NCBI accession ID for the top hit')	# catch URL connection errors
				return mapTuple

		# chromParse: uses Entrez to parse DNA sequence out of chromosome according to BLAST output coordinates
		# input: Chromosome GI number (0), Start coordinate (1), End coordinate (2)
		# output: DNA sequence (0)
		# NOTE: Biopython BLAST documentation: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc81
		def chromParse(GI,start,end):
			handle = Entrez.efetch(db="nucleotide", 
						id=GI, 
						rettype="fasta", 
						strand=1, 
						seq_start=start, 
						seq_stop=end)
			record = SeqIO.read(handle, "fasta")
			Seq = record.seq
			dna_seq = Seq.tostring()
			handle.close()
			return dna_seq

		# tryNCBI: tries does uniprot --> NBCI nucleotide ID fetch
		# input: uniprot ID
		# output: tuple (uniprot ID, NCBI nucleotid ID,DNA ORF matches to standard output
		# NOTE: Programmatic Uniprot Access: http://www.uniprot.org/faq/28#id_mapping_examples
		def tryNCBI(uniprotID):

			url = 'http://www.uniprot.org/mapping/'			# Base URL
			# DB Identifiers: http://www.uniprot.org/faq/28#id_mapping_examples
			mapTuple = ()
			params = {
			'from':'ACC',									# UniProtKB AC, direction "to"	
			'to':'REFSEQ_NT_ID',							# RefSeq Nucleotide, direction "both"
			'format':'tab',									# tabular output format
			'query':uniprotID								# protein query uniprot IDs
			}
			try:
				data = urllib.urlencode(params)				# simultaneous URL queries
				request = urllib2.Request(url, data)
				response = urllib2.urlopen(request)			
				ncbiIDs = response.read(200000)
				ncbiIDs = ncbiIDs.split()					# Get out ids
				if len(ncbiIDs) > 2:						# If it is not empty, proceed
					ncbiIDs = ncbiIDs[2:]
					nucID = ncbiIDs[1]
					sec = nucID[1:2]
					if sec == 'C' or sec == 'G' or sec == 'T' or sec == 'W' or sec == 'Z' or sec == 'S':
						return 'pass'
					handle = Entrez.efetch(db="nucleotide", id=nucID,rettype="fasta")	# use e-fetch to get sequence
					record = SeqIO.read(handle, "fasta")		# use SeqIO to get out the fasta DNA/mRNA sequence
					Seq = record.seq
					mapTuple = (uniprotID,nucID,Seq.tostring(),'NCBI')
					handle.close()
			except:
				# sys.stderr.write('problem reading: '+url+' or performing Entrez efetch')	# catch URL connection errors
				mapTuple = ()
			return mapTuple

		# tryEMBL: tries does uniprot --> EBI nucleotide ID fetch
		# input: uniprot ID
		# output: tuple (uniprot ID, NCBI nucleotide ID,DNA ORF matches to standard output
		# NOTE: Programmatic Uniprot Access: http://www.uniprot.org/faq/28#id_mapping_examples
		def tryEMBL(uniprotID):
			
			url = 'http://www.uniprot.org/mapping/'			# Base URL
			# DB Identifiers: http://www.uniprot.org/faq/28#id_mapping_examples
			mapTuple = ()
			params = {
			'from':'ACC',									# UniProtKB AC, direction "to"	
			'to':'EMBL',									# EBI nucleotide, direction "both"
			'format':'tab',									# tabular output format
			'query':uniprotID								# protein query uniprot IDs
			}
			try:
				data = urllib.urlencode(params)				# simultaneous URL queries
				request = urllib2.Request(url, data)
				response = urllib2.urlopen(request)			
				ebiIDs = response.read(200000)
				# print ebiIDs
				ebiIDs = ebiIDs.split()					# Get out ids
				response.close()
				if len(ebiIDs) > 2:						# If it is not empty, proceed
					ebiIDs = ebiIDs[2:]
					nucID = ebiIDs[1]
				ebiURL = 'http://www.ebi.ac.uk/ena/data/view/'+nucID+'&display=fasta'
				emblFetch = urllib2.urlopen(ebiURL)
				emblFasta = emblFetch.read()
				emblLines = emblFasta.split('\n')
				emblFetch.close()
				emblLines.pop(0)						# Remove fasta description > ...
				DNA = ''.join(emblLines)				
				DNA = DNA.strip()						# Remove whitespace
				mapTuple = (uniprotID,nucID,DNA,'EMBL')
			
			except:
				# sys.stderr.write('problem reading: '+url)	# catch URL connection errors
				mapTuple = ()

			return mapTuple
				
		# getTaxon: get common name for organism corresponding to uniprot ID
		# input: uniprot ID
		# output: tuple
		# NOTE: Programmatic Uniprot Access: http://www.uniprot.org/faq/28#id_mapping_examples
		def getTaxon(uniprotID):
			try:
				page = urllib.urlopen('http://www.uniprot.org/uniprot/'+uniprotID+'.fasta').read()
				pageLines = page.split('\n')
				carrot = pageLines.pop(0)						# Remove fasta description > ...
				match = re.search(r'OS=(\w+\s+\w+)\s',carrot)
				taxName = match.group(1)
				return taxName
			except:			
				# print >>fh,'problem with Uniprot URL connection, or regular expression parsing.' 	# catch URL connection errors
				# sys.stderr.write('problem with Uniprot URL connection, or regular expression parsing.')
				return 'null'

		table = 1										# Translation Table for Biopython, see above
		nucID = 'null'									# Set all default values to null for error catching
		nucSequence = 'null'
		aligned = 'null'
		pID = 'null'
		startORF = 'null'
		startORF = 'null'
		endORF = 'null'
		source = 'null'
		mappedDNA = 'null'
		proSeq = 'null'
		doEMBL = True
		quadTuple = tryNCBI(uniprotID)
		if quadTuple == 'pass':
			doEMBL = False 		# Skip directly to TBLASTN because you need the GI and genomic coordinates
		elif len(quadTuple) == 0:
			quadTuple = tryEMBL(uniprotID)
		if len(quadTuple) != 0 and doEMBL:
			proTuples = orfLength(uniprotID)				# proTuples will tell you the length of the protein 
			min_pro_len = proTuples[0]						# protein length
			proSeq = proTuples[1]							# protein sequence
			nucID = quadTuple[1]
			nucSequence = quadTuple[2]
			dna_seq = Seq(nucSequence)
			source = quadTuple[3]
			(startORF, endORF, aligned, pID, mappedDNA) = orfFinder(dna_seq, table, min_pro_len, proSeq,source)	# Subroutine calls to determine
																			# if matching ORF is present in the
																		# nucleotide sequence and find the
																			# start and stop indices.
		else: 
			taxon = getTaxon(uniprotID)
			if taxon != 'null':
				proTuples = orfLength(uniprotID)				# proTuples will tell you the length of the protein 
				min_pro_len = proTuples[0]						# protein length
				proSeq = proTuples[1]							# protein sequence
				quadTuple = taxBLASTn(uniprotID,taxon,proSeq)
				if len(quadTuple) != 0:
					nucID = quadTuple[1]
					nucSequence = quadTuple[2]
					dna_seq = Seq(nucSequence)
					source = quadTuple[3]
					aligned = quadTuple[4]						# Get aligned from BLAST rather than AlignIO
					(startORF, endORF, dummy, pID, mappedDNA) = orfFinder(dna_seq, table, min_pro_len, proSeq,source)
																				# if matching ORF is present in the
																				# nucleotide sequence and find the
																				# start and stop indices.
		outputTuple = (uniprotID, nucID, source, proSeq, nucSequence, aligned, pID, startORF, endORF, mappedDNA)
		return outputTuple
			
if __name__ == '__main__':
	sys.exit(1)