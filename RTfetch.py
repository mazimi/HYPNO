#!/usr/bin/python -tt
# Copyright Nima Emami, 2012

from Bio import Entrez, SeqIO, AlignIO
from xml.etree.ElementTree import ElementTree
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
import urllib,urllib2, sys, re, xml

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
	def getSeqID(self, uniprotID):
		
		# orfFinder: 6 frame ORF search on DNA/mRNA sequence (Adapted & Modified from BioPython)
		# input: Nucleic Acid Sequence (0), Translation Table Number (1), ORF lenght lower bound (2),
		#			Protein Sequence for alignment confirmation (3)
		# output: Start Index of Hit (0), End Index of Hit (1), flag = 1 (2) OR flag = 0 (0)
		# NOTE: Source code: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc224
			# NOTE: Translation tables: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
		# May need to be modified based on taxonomic classifications, currently set to classical genetic
		# code, although it does not seem to be using this (picking up alternative start codons...)
		def orfFinder(seq, trans_table, min_protein_length, proSeq, source):
			indices = []
			flag = 0														# If set true, match was found
			bestSoFar = 0
			startSoFar = 0
			endSoFar = 0
			pidSoFar = 0
			alignedSoFar = ''
			seq_len = len(seq)
			for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:	# Forward vs. Reverse frames
				for frame in range(3):											# 3 reading frames for each
					trans = str(nuc[frame:].translate(trans_table))			# Translation based on trans_table
					trans_len = len(trans)
					aa_start = 0
					aa_end = 0
					while aa_start < trans_len:
						aa_end = trans.find("*", aa_start)					# Find Stop codon
						if aa_end == -1:
							aa_end = trans_len
						if aa_end-aa_start == min_protein_length:	# If this is the same number of AA's as the Protein
							(percentID, aligned) = getpID(trans[aa_start:aa_end].lower(),proSeq.lower())
							if int(percentID) > 90:
								if strand == 1:
									start = frame+aa_start*3
									end = min(seq_len,frame+aa_end*3+3)
								else:
									start = seq_len-frame-aa_end*3-3
									end = seq_len-frame-aa_start*3
								return (start,end,aligned,int(percentID))
							elif int(percentID) > bestSoFar:	# If there is greater than 97 % ID
												# (See isoform discussion: http://www.uniprot.org/faq/30)
								bestSoFar = int(percentID)

								if strand == 1:
									start = frame+aa_start*3
									end = min(seq_len,frame+aa_end*3+3)
								else:
									start = seq_len-frame-aa_end*3-3
									end = seq_len-frame-aa_start*3                        
								# indices.append(start, end, strand,# trans[aa_start:aa_end]))
								startSoFar = start
								endSoFar = end
								alignedSoFar = aligned
								pidSoFar = percentID
						aa_start = aa_start+1
			return (startSoFar,endSoFar,alignedSoFar,pidSoFar)

		# pID: pairwise alignment subroutine (Adapted and modified from below source)
		# input: 2 nucleotide sequences (0,1)
		# output: integer value between 0 and 100 (0)
		# NOTE: Source code: http://biostumblematic.wordpress.com/2009/06/15/measuring-identities-of-aligned-protein-sequences-with-biopython/
		# NOTE: Prints to temporary file handle, may want to find a way to avoid this if possible
		def getpID(seq1,seq2):											
			fasta_temp = open('temp.fasta', 'w')					
			print >>fasta_temp, '> seq1\n'+seq1+'\n> seq2\n'+seq2
			fasta_temp.close()
			fastas = open('temp.fasta', 'rU')
			alignment = AlignIO.read(fastas, "fasta")			# AlignIO BioPerl Module
		 	aligned = ''
			j=0 # counts positions in first sequence
			i=0 # counts identity hits
			for record in alignment[1:]:
				aligned = record.seq
				for amino_acid in record.seq:
					if amino_acid == '-':
						pass
					else:
						if amino_acid == alignment[0].seq[j]:
							i += 1
					j += 1
				j = 0
			seq = str(record.seq)
			gap_strip = seq.replace('-', '')
			percent = 100*i/len(gap_strip)
			return (str(percent),str(aligned).upper())

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
				# print >>fh,'Q: '+qseq
				# print >>fh,"M: "+midseq
				# print >>fh,"H: "+hseq
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
					handle = Entrez.efetch(db="nucleotide", id=nucID,rettype="fasta")	# use e-fetch to get sequence
					record = SeqIO.read(handle, "fasta")		# use SeqIO to get out the fasta DNA/mRNA sequence
					Seq = record.seq
					mapTuple = (uniprotID,nucID,Seq.tostring(),'NCBI')
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
		nucID = 'null'
		nucSequence = 'null'
		aligned = 'null'
		pID = 'null'
		startORF = 'null'
		startORF = 'null'
		endORF = 'null'
		source = 'null'
		quadTuple = tryNCBI(uniprotID)
		if len(quadTuple) == 0:
			quadTuple = tryEMBL(uniprotID)
		if len(quadTuple) != 0:
			proTuples = orfLength(uniprotID)				# proTuples will tell you the length of the protein 
			min_pro_len = proTuples[0]						# protein length
			proSeq = proTuples[1]							# protein sequence
			nucID = quadTuple[1]
			nucSequence = quadTuple[2]
			dna_seq = Seq(nucSequence)
			source = quadTuple[3]
			(startORF, endORF, aligned, pID) = orfFinder(dna_seq, table, min_pro_len, proSeq,source)	# Subroutine calls to determine
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
					(startORF, endORF, dummy, pID) = orfFinder(dna_seq, table, min_pro_len, proSeq,source)
																				# if matching ORF is present in the
																				# nucleotide sequence and find the
																				# start and stop indices.
		outputTuple = (uniprotID, nucID, source, proSeq, nucSequence, aligned, pID, startORF, endORF)
		return outputTuple
			
if __name__ == '__main__':
	sys.exit(1)