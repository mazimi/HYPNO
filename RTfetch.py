#!/usr/bin/python -tt
# Copyright Nima Emami, 2012

from Bio import Entrez, SeqIO, AlignIO
from xml.etree.ElementTree import ElementTree
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
import urllib,urllib2, sys, re, xml, os

# Installation Requirements: ElementTree module, EMBOSS package, and BioPython for functionality

class RTfetch:

	# __init__: Constructor method, delegates nucleotide sequence fetch
	# input: self
	# output: none
	def __init__(self):
		Entrez.email = 'timanoed154@yahoo.com'				# mandatory NCBI email for reference

	# getSeqID: Delegates nucleotide sequence fetch
		# input: self (0), uniprotID (1)
		# output: 10 element tuple of nucleotide fetch details (see below)
		# NOTE: Programmatic Uniprot Access: http://www.uniprot.org/faq/28#id_mapping_examples
			# Tuple indices:
			# (0) Uniprot ID query
			# (1) Nucleotide ID (default = 'null')
			# (2) Database source (EMBL, NCBI, TBLASTN, default = 'null')
			# (3) Protein sequence (default = 'null')
			# (4) DNA sequence (default = 'null')
			# (5) Aligned protein hit (default = 'null')
			# (6) Percent identity (default = 'null')
			# (7) ORF starting index (default = 'null')
			# (8) ORF ending index (default = 'null')
			# (9) Mapped DNA sequence (default = 'null')
	def getSeqID(self, uniprotID):

		# orfFinder: 6 frame ORF search on nucleotide sequence (adapted & modified from BioPython Cookbook, like below)
		# input: nucleic acid sequence (0), translation table integer value (1), target protein length (2),
		#			protein sequence for pairwise alignment comparisons (3)
		# output: Start Index of Hit (0), End Index of Hit (1), flag = 1 (2) OR flag = 0 (0)
		# NOTE: adapted based on BioPython code: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc224
		def orfFinder(seq, trans_table, targetLength, proSeq):
			bestSoFar = 0 					# Set defaults
			gapsSoFar = targetLength*3		# ^
			startSoFar = 0 					# |
			pStartSoFar = 0 				# |
			pEndSoFar = 0 					# |
			endSoFar = 0 					# |
			pidSoFar = 0 					# |
			alignedSoFar = '' 				# |
			querySoFar = ''					# |
			mappedSeq = '' 					# |
			seq_len = len(seq)				# length of original protein sequence, for discriminatory ORF finding
			for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:	# Forward vs. Reverse frames
				for frame in range(3):											# 3 reading frames for each
					trans = str(nuc[frame:].translate(trans_table))
					trans_len = len(trans)
					aa_start = 0
					aa_end = 0
					while aa_start < trans_len:
						aa_end = trans.find("*", aa_start)	# find ter (stop) codon, labeled as '*'
						if aa_end == -1:
							aa_end = trans_len
						# ** current heuristic: only look for nucleotide sequences which code for proteins 3/4ths to 5/4ths of target size
						# -->  May want to refine this at some point
						if 1.25*targetLength >= aa_end-aa_start >= 0.75*targetLength:
							(percentID, aligned, alignedQuery, gaps) = getpID(trans[aa_start:aa_end].lower(),proSeq.lower())
							if gaps == 0:				# If it is a perfect match
								if strand == 1:						# Forward frame
									start = max(0,frame+aa_start*3)			# start index (w/ max because aa_start + frame may be less than 0)
									end = min(seq_len,frame+aa_end*3+3)		# end index (w/ +3 at end because frame goes negative instead of positive)
									mappedSeq = mapDNA(start,end,seq,aligned,alignedQuery)	# Call DNA mapping subroutine to map DNA to aligned protein to query
								else:								# Reverse frame
									start = max(0,seq_len-frame-aa_end*3-3)	# Backwards case, so everything is effectively reversed
									end = seq_len-frame-aa_start*3 			# end index is now at the beginning, so subtract off AA start
									mappedDNA = mapDNA(start,end,seq,aligned)	# Call DNA mapping subroutine to map DNA to aligned protein to query
									mappedSeq = Seq(mappedDNA)					# RevComp the DNA sequence
									mappedSeq = mappedSeq.reverse_complement()
								return (end,start,aligned,int(percentID),str(mappedSeq)) # Flip Start and End indices because of start and stop index
							elif gaps < gapsSoFar:	# Trying to maximize percent identity
													# (See isoform discussion: http://www.uniprot.org/faq/30)
								gapsSoFar = gaps
								if strand == 1:
									start = max(0,frame+aa_start*3)
									end = min(seq_len,frame+aa_end*3+3)
								else:
									start = max(0,seq_len-frame-aa_end*3-3)
									end = seq_len-frame-aa_start*3                     
								startSoFar = start
								endSoFar = end
								alignedSoFar = aligned
								querySoFar = alignedQuery
								pidSoFar = percentID
						aa_start = aa_start+1
			mappedDNA = mapDNA(startSoFar,endSoFar,seq,alignedSoFar,querySoFar)
			return (startSoFar,endSoFar,alignedSoFar,int(pidSoFar),mappedDNA)

		# mapDNA: mapping DNA sequence to the Needleman-Wunsch aligned protein to the original protein query
		# input: start index (0), end index (1), dna (2), aligned translated protein (3), aligned original query (4)
		# output: mapped DNA sequence
		def mapDNA(startIndex,endIndex,dna,alignedProtein,alignedQuery):
			dna = str(dna)
			mapped = ''
			dnaCounter = startIndex
			protCounter = 0
			print "print"
			Start alignedProtein
			print
			print alignedQuery
			print
			while(protCounter < len(alignedProtein)):	# While you have not surpassed the length of the aligned protein sequence
				if alignedProtein[protCounter] == '-':	# If it is a gap in the translated sequence, then add '***' but
					mapped = mapped + '***'					# don't skip forward in the DNA sequence yet
				elif alignedQuery[protCounter] == '-':	# If it is a gap in the original sequence, this corresponds to an insertion
					pass									# so don't insert any gaps
				elif alignedProtein[protCounter] != alignedQuery[protCounter]:
					mapped = mapped + '***'				# Insert three gaps '***' in the DNA alignment and also skip ahead in the DNA
					dnaCounter = dnaCounter + 3 			# because the current AA was mismatched
				else:
					mapped = mapped + dna[dnaCounter:dnaCounter+3]	# Otherwise, insert the next three nucleic acids
					dnaCounter = dnaCounter + 3 					# and increment both counters
				protCounter = protCounter + 1 			# In all cases, move foward in the aligned protein sequences
			return mapped

		# pID: pairwise alignment subroutine (Adapted and modified from below source)
		# input: 2 nucleotide sequences (0,1)
		# output: integer value between 0 and 100 (0), aligned translated protein sequence (1)
		# NOTE: Uses EMBOSS package needle executable call
		# NOTE: Creates, prints to, and reads from local file tempAlign.needle
		def getpID(target,template):											
			seqA = open('tempA.fasta', 'w')			# Print sequences to temporary files for EMBOSS needle call
			seqB = open('tempB.fasta', 'w')				
			print >>seqA, '> seqA\n'+target+'\n'
			print >>seqB, '> seqB\n'+template+'\n'
			seqA.close()							# Close temporary file handle references
			seqB.close()
			os.system('needle -asequence tempA.fasta -sprotein1 -bsequence tempB.fasta -sprotein2 -gapopen 10 -gapextend 0.5 -outfile tempAlign.needle -auto')
			needle = open('tempAlign.needle','rU')
			alignment = AlignIO.read(needle,"emboss")# AlignIO BioPython Module reads out EMBOSS globally aligned sequences
			i=0 									 # Global alignment --> only 1 counter necessary for both sequences
			counter = 0 							 # Global counter
			gaps = 0
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
				if topAA != bottomAA:
					gaps = gaps + 1
					pass
				else:
					i = i + 1
				counter = counter + 1
			percent = 100*i/len(seq0)
			os.remove('tempA.fasta')	# Remove these temporary files that you no longer need
			os.remove('tempB.fasta')
			return (str(percent),seq0.upper(),seq1.upper(), gaps)

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
		# output: DNA sequence corresponding to that NCBI ID
		def getSeq(NCBI_ID):
			try:
				handle = Entrez.efetch(db="nucleotide", id=NCBI_ID,rettype="fasta")	# use e-fetch to get sequence
				record = SeqIO.read(handle, "fasta")		# use SeqIO to get out the fasta DNA/mRNA sequence
				Seq = record.seq
				dna_seq = Seq.tostring()
			except:
				# print >>fh,'problem reading: '+url+' or performing Entrez efetch' # catch URL connection errors
				# sys.stderr.write('problem reading: '+url+' or performing Entrez efetch')	# catch URL connection errors
				dna_seq = 'null'
			return dna_seq

		# taxBLASTn: tries does uniprot --> NBCI nucleotide ID fetch
		# input: uniprot ID (0), common taxonomic organism name (1)
		# output: tuple (uniprot ID, NCBI nucleotide ID,DNA ORF matches to standard output
		# NOTE: BioPython BLAST documentation: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc81
		def taxBLASTn(uniprotID,taxName,proSeq):
			mapTuple = ()
			try:
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
				prefix = accessionNCBI[0:2]
				# ** NOTE: based on following keys:
				# 1) http://www.ncbi.nlm.nih.gov/RefSeq/key.html
				# 2) http://www.ncbi.nlm.nih.gov/Sequin/acc.html
				# TODO: improve robustness
				flag = prefix == 'AC' or prefix == 'NC' or prefix == 'NG' or prefix == 'NT' or prefix == 'NW'\
							or prefix == 'NZ' or prefix == 'NS' or prefix == 'AP' or prefix == 'BS' or prefix == 'AL'\
							or prefix == 'BX' or prefix == 'CR' or prefix == 'CT' or prefix == 'CU' or prefix == 'FP'\
							or prefix == 'FQ' or prefix == 'FR' or prefix == 'AE' or prefix == 'CP' or prefix == 'CY'\
							or prefix == 'AM' 
				if flag:
					dna_seq = chromParse(GI,Hit_from,Hit_to)
				else:
					dna_seq = getSeq(accessionNCBI)
				if dna_seq != '':
					mapTuple = (uniprotID,accessionNCBI,dna_seq,'TBLASTN',midseq)
					return mapTuple
				else:
					return mapTuple
			except:
				return mapTuple

		# chromParse: uses Entrez to parse DNA sequence out of chromosome according to BLAST output coordinates
		# input: chromosome GI number (0), start coordinate (1), end coordinate (2)
		# output: DNA sequence (0)
		# NOTE: BioPython BLAST documentation: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc81
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

		# tryNCBI: tries uniprot --> NBCI nucleotide ID fetch
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
					prefix = nucID[0:2]
					# ** NOTE: based on following keys:
					# 1) http://www.ncbi.nlm.nih.gov/RefSeq/key.html
					# 2) http://www.ncbi.nlm.nih.gov/Sequin/acc.html
					# TODO: improve robustness
					if prefix == 'AC' or prefix == 'NC' or prefix == 'NG' or prefix == 'NT' or prefix == 'NW'\
							or prefix == 'NZ' or prefix == 'NS' or prefix == 'AP' or prefix == 'BS' or prefix == 'AL'\
							or prefix == 'BX' or prefix == 'CR' or prefix == 'CT' or prefix == 'CU' or prefix == 'FP'\
							or prefix == 'FQ' or prefix == 'FR' or prefix == 'AE' or prefix == 'CP' or prefix == 'CY'\
							or prefix == 'AM':
							return 'pass'
					handle = Entrez.efetch(db="nucleotide", id=nucID,rettype="fasta")	# use e-fetch to get sequence
					record = SeqIO.read(handle, "fasta")		# use SeqIO to get out the fasta DNA/mRNA sequence
					Seq = record.seq
					mapTuple = (uniprotID,nucID,Seq.tostring(),'NCBI')
					handle.close()
			except:
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
				mapTuple = ()

			return mapTuple

		# getTaxon: get common name for organism corresponding to uniprot ID
		# input: uniprot ID
		# output: organism common name
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
				return 'null'

		# getGeneticCode: for deducing organismically specific nucleotide translation table based on UniProt taxonomic lineage of query
		# input: uniprot ID
		# output: integer for genetic code 
		# NOTE: Translation tables: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
		def getGeneticCode(uniprotID):
			try:
				value = 1 		# Preset: default genetic code integer value
				TaxonomicLineage = []
				ns = '{http://uniprot.org/uniprot}'
				pageXML = urllib.urlopen('http://www.uniprot.org/uniprot/'+uniprotID+'.xml').read()
				tree = xml.etree.ElementTree.fromstring(pageXML)
				items = tree.getiterator(ns+'uniprot')
				item = items[0]
				taxa = item.find(ns+'entry').find(ns+'organism').find(ns+'lineage').findall(ns+'taxon')
				for taxon in taxa:
					TaxonomicLineage.append(taxon.text)
				for taxon in TaxonomicLineage:
					# NOTE: need to be careful with these, alternative translation may apply to some
					# sequences from the specific taxa but not others (i.e. only MITOCHONDRIAL sequences...)
					# if taxon == 'Vertebrata':
					# 	value = 2
					# elif taxon == 'Saccharomyces' or taxon == 'Candida' or taxon == 'Hansenula' or taxon == 'Kluyveromyces':
					# 	value = 3
					# elif taxon == 'Entomoplasmatales' or taxon == 'Mycoplasmatales':
					# 	value = 4
					# elif taxon == 'Nematoda' or taxon == 'Mollusca' or taxon == 'Arthropoda':
					# 	value = 5
					# elif taxon == 'Ciliata' or taxon == 'Dasycladaceae' or taxon == 'Diplomonadida':
					# 	value = 6
					if taxon == 'Asterozoa' or taxon == 'Echinozoa' or taxon ==  'Rhabditophora':
						value = 9
					elif taxon == 'Bacteria' or taxon == 'Archaea':
						value = 11
					# elif taxon == 'Ascidia':
					# 	value = 13
					# elif taxon == 'Chlorophyceae':
					# 	value = 16
					# elif taxon == 'Scenedesmus':
					# 	value = 22
					# elif taxon == 'Trematoda':
					# 	value = 21
					# elif taxon == 'Thraustochytrium':
					# 	value = 23
					# elif taxon == 'Pterobranchia':
					# 	value = 24
				return value
			except:			
				return 1 		# Exceptional Case: throw default genetic code integer value

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
		table = getGeneticCode(uniprotID) 				# Translation Table for Biopython, default = 1 (standard genetic code)
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
			(startORF, endORF, aligned, pID, mappedDNA) = orfFinder(dna_seq, table, min_pro_len, proSeq)
			outputTuple = (uniprotID, nucID, source, proSeq, nucSequence, aligned, pID, startORF, endORF, mappedDNA)
			return outputTuple
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
					(startORF, endORF, dummy, pID, mappedDNA) = orfFinder(dna_seq, table, min_pro_len, proSeq)
			outputTuple = (uniprotID, nucID, source, proSeq, nucSequence, aligned, pID, startORF, endORF, mappedDNA)
			return outputTuple

if __name__ == '__main__':
	sys.exit(1)