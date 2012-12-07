#!/usr/bin/python -tt
# Copyright Nima Emami, 2012

from Bio import Entrez, SeqIO, AlignIO
from xml.etree.ElementTree import ElementTree
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
import os, re, subprocess, sys, xml, urllib,urllib2
from Bio.Emboss.Applications import NeedleCommandline
from HTMLParser import HTMLParser

# Dependencies: EMBOSS needle, Biopython, internet connection

class RTfetch:

	# __init__: Constructor method, delegates nucleotide sequence fetch
	# input: self
	# output: none
	def __init__(self):
		# mandatory NCBI Entrez email for tBLASTn reference
		Entrez.email = 'bpgHYPNO@gmail.com'

	# getSeqID: Delegates nucleotide sequence fetch
		# input: self (0), UniProt accession (1)
		# output: 11 element tuple of nucleotide fetch details (see below)
		# NOTE: Programmatic Uniprot Access: http://www.uniprot.org/faq/28#id_mapping_examples
			# Tuple indices (all defaults set to 'null', except for index 0):
			# (0) UniProt ID query
			# (1) Nucleotide sequence accession
			# (2) Database source (EMBL, NCBI, TBLASTN)
			# (3) UniProt amino acid sequence
			# (4) Nucleic acid sequence
			# (5) Predicted protein sequence, with gaps
			# (6) Predicted protein % identity to index (3) sequence
			# (7) ORF starting index
			# (8) ORF ending index
			# (9) Mapped DNA sequence, with gaps
			# (10) Error message
			# (11) Debug Info
	def getSeqID(self, uniprotID):

		# orfFinder: nucleotide ORF search in all 6 frames
		# input: nucleic acid sequence (0), translation table integer value (1), target protein length (2),
		#			protein sequence for pairwise alignment comparisons (3)
		# output: ORF start index (0), ORF end index (1), aligned predicted protein sequence (2), predicted protein
		# 			pID relative to UniProt sequence (3), sequence of matching nucleotides mapped onto UniProt sequence (4)
		# NOTE: adapted based on BioPython code: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc224
		def orfFinder(seq, trans_table, targetLength, proSeq):
			bestSoFar, gapsSoFar, startSoFar, pStartSoFar, pEndSoFar, endSoFar, pidSoFar = 0, targetLength*3, 0, 0, 0, 0, 0
			alignedSoFar, querySoFar, mappedSeq, seq_len = '', '', '', len(seq)
			for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:	# Forward vs. Reverse frames
				for frame in range(3):										# 3 reading frames for each
					trans = str(nuc[frame:].translate(trans_table))
					trans_len, aa_start, aa_end = len(trans), 0, 0
					while aa_start < trans_len:
						aa_end = trans.find("*", aa_start)					# Find stop codon ('*' character)
						if aa_end == -1:
							aa_end = trans_len
						# Heuristic: look for ORFs which code for proteins 1/2x to 3/2x the target size
						if 1.50*targetLength >= aa_end-aa_start >= 0.50*targetLength:
							(percentID, aligned, alignedQuery, gaps) = getpID(trans[aa_start:aa_end].lower(),proSeq.lower())
							if gaps == 0:									# Perfect match
								if strand:									# Forward frame
									# start index (max() because aa_start + frame may be less than 0)
									start = max(0,frame+aa_start*3)
									# end index (+3 because frame goes negative instead of positive)
									end = min(seq_len,frame+aa_end*3+3)
									# "align" ORF DNA to aligned predicted protein		
									mappedSeq = mapDNA(start,end,seq,aligned,alignedQuery)	
								else:										# Reverse frame
									# Backwards case, so everything is effectively reversed
									start = max(0,seq_len-frame-aa_end*3-3)
									# end index is now at the beginning, so subtract off aa_start
									end = seq_len-frame-aa_start*3 			
									mappedDNA = mapDNA(start,end,seq,aligned)	
									mappedSeq = Seq(mappedDNA).reverse_complement()
								# Flip start and end indices because of start/stop index
								return (end,start,aligned,percentID,str(mappedSeq))
							elif gaps < gapsSoFar:	# Trying to maximize percent identity
								# (See isoform discussion: http://www.uniprot.org/faq/30)
								gapsSoFar = gaps
								if strand == 1:
									start, end = max(0,frame+aa_start*3), min(seq_len,frame+aa_end*3+3)
								else:
									start, end  = max(0,seq_len-frame-aa_end*3-3), seq_len-frame-aa_start*3
								# Update information for best hit seen at this iteration
								startSoFar, endSoFar, alignedSoFar, querySoFar, pidSoFar = start, end, aligned, alignedQuery, percentID
						aa_start += 1
			mappedDNA = mapDNA(startSoFar,endSoFar,seq,alignedSoFar,querySoFar)
			return startSoFar,endSoFar,alignedSoFar,pidSoFar,mappedDNA

		# mapDNA: mapping DNA sequence to the Needleman-Wunsch aligned predicted protein, relative to UniProt sequence
		# input: start index (0), end index (1), DNA sequence (2), aligned translated protein (3), aligned original query (4)
		# output: mapped DNA sequence
		def mapDNA(startIndex,endIndex,dna,alignedProtein,alignedQuery):
			dna, mapped, dnaCounter, protCounter = str(dna), '', startIndex, 0
			while(protCounter < len(alignedProtein)):	# While you have not surpassed the length of the aligned protein sequence
				if alignedProtein[protCounter] == '-':	# If it is a gap in the translated sequence, then add '***' but
					mapped += '***'							# don't skip forward in the DNA sequence yet
				elif alignedQuery[protCounter] == '-':	# If it is a gap in the original sequence, this corresponds to an insertion
					pass									# so don't insert any gaps
				elif alignedProtein[protCounter].upper() != alignedQuery[protCounter].upper():
					mapped += '***'				# Insert three gaps '***' in the DNA alignment and also skip ahead in the DNA
					dnaCounter += 3 			# because the current AA was mismatched
				else:
					mapped = mapped + dna[dnaCounter:dnaCounter+3]	# Otherwise, insert the next three nucleic acids
					dnaCounter += 3 					# and increment both counters
				protCounter += 1 						# In all cases, move foward in the aligned protein sequences
			return mapped

		def getpID(target,template):
			# Write sequences to temporary files for EMBOSS needle call											
			with open('tempA.fasta', 'w') as seqA:			
				with open('tempB.fasta', 'w') as seqB:			
					seqA.write('> seqA\n'+target+'\n')
					seqB.write('> seqB\n'+template+'\n')
			os.system('needle -asequence tempA.fasta -sprotein1 -bsequence tempB.fasta -sprotein2 -gapopen 10 -gapextend 0.5 \
							-outfile tempAlign.needle -auto')
			with open('tempAlign.needle','rU') as needle:
				alignment = AlignIO.read(needle,"emboss")
				i, counter, gaps = float(0), 0, 0
				sequence0, sequence1= alignment[0], alignment[1]
				seq0, seq1 = str(sequence0.seq), str(sequence1.seq)
				chars0, chars1 = list(seq0), list(seq1)
				while counter < len(chars0):
					topAA = chars0[counter]
					bottomAA = chars1[counter]
					# Considers gaps and mismatches in both sequences for computing percent identity
					if topAA.upper() != bottomAA.upper():
						gaps += 1
						pass
					else:
						i += float(1)
					counter += 1
				percent = float(100)*i/float(len(seq0))
				# Dispose of temporary FASTA files
			os.remove('tempA.fasta')
			os.remove('tempB.fasta')
			os.remove('tempAlign.needle')
			return float(percent),seq0.upper(),seq1.upper(), gaps

		# orfLength: determines length of ORFs based on protein sequence length via uniprot query
		# input: uniprot ID (0)
		# output: tuple of (length of protein, protein sequence)
		def orfLength(protID):
			page = urllib.urlopen('http://www.uniprot.org/uniprot/'+protID+'.fasta').read()
			pageLines = page.split('\n')
			pageLines.pop(0)
			merged = ''.join(pageLines).strip()				
			return len(merged),merged

		# getSeq: retrieve nucleotide sequence from Entrez using EMBL / NCBI accession
		# input: accession (0)
		# output: DNA sequence (0)
		def getSeq(accession):
			try:
				handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta")
				record = SeqIO.read(handle, "fasta")
				dnaSeq = record.seq.tostring()
				handle.close()
			except:
				return 'null'
			return dnaSeq

		# taxBLASTn: taxonomically specific tBLASTn maps uniprot --> NBCI nucleotide ID
		# input: UniProt ID (0), Scientific name (1)
		# output: tuple (UniProt accession, NCBI nucleotide accession, DNA ORF matches)
		# NOTE: BioPython BLAST documentation: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc81
		def taxBLASTn(uniprotID,taxName,proSeq):
			mapTuple, debug = (), ''
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
				debug += 'tBLASTn hit accession and match indices: '+str(accessionNCBI)+' ('+str(Hit_to)+', '+str(Hit_from)+')\n'
				dna_seq = chromParse(GI,Hit_from,Hit_to)
				if dna_seq != '':
					mapTuple = (accessionNCBI,dna_seq,'TBLASTN',midseq)
					return mapTuple, debug
				else:
					return mapTuple, debug
			except:
				return mapTuple, debug

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
			dna_seq = record.seq.tostring()
			handle.close()
			return dna_seq

		def accessionHTTP(to, accession):
			url, nucID = 'http://www.uniprot.org/mapping/', ''			# Base URL
			# DB Identifiers: http://www.uniprot.org/faq/28#id_mapping_examples
			params = {
				'from':'ACC',								# UniProtKB AC, direction "to"	
				'to':to,									# ID source
				'format':'tab',								# tabular output format
				'query':accession							# protein query uniprot IDs
			}
			try:
				data = urllib.urlencode(params)				# simultaneous URL queries
				request = urllib2.Request(url, data)
				response = urllib2.urlopen(request)			
				IDs = response.read(200)
				IDs = IDs.split()					# Get out ids
				response.close()
				if len(IDs) > 2:					# If it is not empty, proceed
					IDs = IDs[2:]
					nucID = IDs[1]
			except:
				pass
			return nucID

		# tryNCBI: tries uniprot --> NBCI nucleotide ID fetch
		# input: uniprot ID
		# output: tuple (uniprot ID, NCBI nucleotid ID,DNA sequence, source)
		# NOTE: Programmatic Uniprot Access: http://www.uniprot.org/faq/28#id_mapping_examples
		def tryNCBI(uniprotID):
			url, nucID, debug = 'http://www.uniprot.org/mapping/', '', ''	# Base URL
			# DB Identifiers: http://www.uniprot.org/faq/28#id_mapping_examples
			mapTuple = ()
			nucID = accessionHTTP('REFSEQ_NT_ID', uniprotID)
			if nucID == '':
				return mapTuple, debug
			try:
				genomic = Entrez.efetch(db="nucleotide", id=nucID, rettype="gb")
				gbHeader = genomic.read(200)
				match = re.search('([\d]+) bp',gbHeader)
				size = int(match.group(1))
				debug += 'For UniProt accession '+uniprotID+', the retrieved NCBI accession is '+nucID+'\n'
				if size > 10000:
						# if genomic --> return 'pass' and move to tBLASTn
						debug += '\tSwitching to tBLASTn mode\n'
						return 'pass', debug
				mapTuple = (nucID,getSeq(nucID),'NCBI')
			except:
				pass
			return mapTuple, debug

		# tryEMBL: tries does uniprot --> EMBL nucleotide ID fetch
		# input: uniprot ID
		# output: tuple (uniprot ID, EMBL nucleotide ID, DNA sequence, source)
		# NOTE: Programmatic Uniprot Access: http://www.uniprot.org/faq/28#id_mapping_examples
		def tryEMBL(uniprotID):
			url, debug = 'http://www.uniprot.org/mapping/', ''			# Base URL
			# DB Identifiers: http://www.uniprot.org/faq/28#id_mapping_examples
			mapTuple = ()
			nucID = accessionHTTP('EMBL_ID', uniprotID)
			if nucID == '':
				debug += 'EMBL accession could not be retrieved, switching to tBLASTn mode.\n'
				return mapTuple, debug
			try:
				genomic = urllib2.Request('http://www.ebi.ac.uk/ena/data/view/'+nucID+'&display=text')
				genText = urllib2.urlopen(genomic)
				IDline = genText.read(200)
				match = re.search('([\d]+) BP',IDline)
				size = int(match.group(1))
				genText.close()
				debug += 'For UniProt accession '+uniprotID+', the retrieved EMBL accession is '+nucID+'\n'
				if size > 10000:
						debug += '\tSwitching to tBLASTn mode\n'
						# if genomic --> return and move to tBLASTn
						return mapTuple, debug
				mapTuple = (nucID,getSeq(nucID),'EMBL')
			except:
				pass
			return mapTuple, debug

		# getTaxon: get common name for organism corresponding to uniprot ID
		# input: uniprot ID
		# output: organism common name
		# NOTE: Programmatic Uniprot Access: http://www.uniprot.org/faq/28#id_mapping_examples
		def getTaxon(uniprotID):
			try:
				page = urllib.urlopen('http://www.uniprot.org/uniprot/'+uniprotID+'.fasta').read()
				pageLines = page.split('\n')
				carrot = pageLines.pop(0)						# Remove fasta description > ...
				match = re.search(r'OS=(\w+\s+[^A-Z]+)\s+[A-Z]',carrot)
				taxName = match.group(1)
				return taxName
			except:			
				return 'null'

		# checkObsolete: check for UniProt accessions that have been removed
		# input: UniProt ID (0)
		# output: tuple containing (obsolete flag, cause of error) (0)
		def checkObsolete(inputID):
			page = urllib.urlopen('http://www.uniprot.org/uniprot/'+inputID)
			read = page.read()
			page.close()
			parser = HTMLUniparse()
			parser.setID(inputID)
			parser.feed(read)
			if parser.obsolete:
				error = 'Obsolete: UniProt Entry '+inputID+"\n"
				if parser.deleted:
					error = error + '\tCause: Entry deleted\n'
				elif parser.demerged:
					error = error + '\tCause: Entry demerged into following UniProt IDs: '+str(parser.parsedIDs)+"\n"
				else:
					error = error + '\tCause: Entry merger\n'
				return (True, error)
			return (False, 'null')

		# getGeneticCode: deducing organismically specific nucleotide translation table from UniProt taxonomic lineage
		# input: uniprot ID (0)
		# output: integer for genetic code (0)
		# NOTE: Translation tables: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
		def getGeneticCode(uniprotID):
			try:
				# Preset: default genetic code integer value
				value, mitochondrial, scientificName = 1, 0, ''
				TaxonomicLineage = []
				ns = '{http://uniprot.org/uniprot}'
				pageXML = urllib.urlopen('http://www.uniprot.org/uniprot/'+uniprotID+'.xml').read()
				tree = xml.etree.ElementTree.fromstring(pageXML)
				item = tree.getiterator(ns+'uniprot')[0]
				taxa = item.find(ns+'entry').find(ns+'organism').find(ns+'lineage').findall(ns+'taxon')
				
				localization = item.find(ns+'entry').find(ns+'geneLocation')
				for key, locale in localization.attrib.iteritems():
					if locale.lower() == 'mitochondrion':
						mitochondrial = True
					
				names = item.find(ns+'entry').find(ns+'organism').findall(ns+'name')
				for name in names:
					for key, nameType in name.attrib.iteritems():
						if nameType.lower() == 'scientific':
							scientificName = name.text
					
				for taxon in taxa:
					TaxonomicLineage.append(taxon.text)
				for taxon in TaxonomicLineage:
					if taxon == 'Vertebrata' and mitochondrial:		
						value = 2
					elif scientificName == 'Saccharomyces cerevisiae' or scientificName == 'Candida glabrata'  \
						or scientificName == 'Hansenula saturnus' or scientificName == 'luyveromyces thermotolerans' \
						and mitochondrial:
						value = 3
					elif taxon == 'Ascaris' or taxon == 'Caenorhabditis' or taxon == 'Bivalvia' \
						or taxon == 'Polyplacophora' or taxon == 'Artemia' or taxon == 'Drosophila' \
						and mitochondrial:
						value = 5
					elif taxon == 'Asterozoa' or taxon == 'Echinozoa' or taxon ==  'Rhabditophora' \
						and mitochondrial:
						value = 9
					elif taxon == 'Euplotidae':
						value = 10
					elif taxon == 'Bacteria' or taxon == 'Archaea':
						value = 11
					elif scientificName == 'Candida albicans' or scientificName == 'Candida cylindracea' \
						or scientificName == 'Candida melibiosica' or scientificName == 'Candida parapsilosis' \
						or scientificName == 'Candida rugosa' and mitochondrial:
						value = 12
					elif taxon == 'Blepharisma':
						value = 15
					elif taxon == 'Chlorophyceae' or scientificName == 'Spizellomyces punctatus' \
						and mitochondrial:
						value = 16
					elif taxon == 'Trematoda' and mitochondrial:
						value = 21
					elif scientificName == 'Scenedesmus obliquus' and mitochondrial:
						value = 22
					elif scientificName == 'Thraustochytrium aureum' and mitochondrial:
						value = 23
					elif taxon == 'Pterobranchia' and mitochondrial:
						value = 24
				return value
			except:			
				return 1 		# Exceptional Case: throw default genetic code integer value

		# Set all default values to null for error catching
		nucID, nucSequence, aligned, pID, startORF, endORF, source, mappedDNA, proSeq, errorMessage  = ('null', ) * 10
		debug = ''
		(isObsolete, error) = checkObsolete(uniprotID)
		if isObsolete:
			errorMessage = error
			outputTuple = (uniprotID, nucID, source, proSeq, nucSequence, aligned, pID, startORF, endORF, mappedDNA, errorMessage, debug)
			return outputTuple
		table = getGeneticCode(uniprotID) 				# Translation Table for Biopython, default = 1 (standard genetic code)
		straight2tBLASTn = False
		triTuple, info = tryNCBI(uniprotID)
		debug += info
		if triTuple == 'pass':
			straight2tBLASTn = True 						# Skip directly to TBLASTN because you need the GI and genomic coordinates
		elif len(triTuple) == 0:							# Otherwise, do EMBL search.
			triTuple, info = tryEMBL(uniprotID)
			debug += info
		if len(triTuple) != 0 and not straight2tBLASTn:
			proTuple = orfLength(uniprotID)				# proTuple will tell you the length of the protein 
			min_pro_len, proSeq  = proTuple
			nucID, nucSequence, source = triTuple
			dna_seq = Seq(nucSequence)
			(startORF, endORF, aligned, pID, mappedDNA) = orfFinder(dna_seq, table, min_pro_len, proSeq)
			outputTuple = (uniprotID, nucID, source, proSeq, nucSequence, aligned, pID, startORF, endORF, mappedDNA, errorMessage, debug)
		else:
			taxon = getTaxon(uniprotID)
			if taxon != 'null':
				proTuple = orfLength(uniprotID)				# proTuple will tell you the length of the protein 
				min_pro_len, proSeq = proTuple					# protein length
				quadTuple, info = taxBLASTn(uniprotID,taxon,proSeq)
				debug += info
				if len(quadTuple) != 0:
					nucID, nucSequence, source, aligned = quadTuple		# Get aligned prot from BLAST rather than AlignIO
					dna_seq = Seq(nucSequence)
					(startORF, endORF, dummy, pID, mappedDNA) = orfFinder(dna_seq, table, min_pro_len, proSeq)
			outputTuple = (uniprotID, nucID, source, proSeq, nucSequence, aligned, pID, startORF, endORF, mappedDNA, errorMessage, debug)
		return outputTuple

# create a subclass and override the handler methods
class HTMLUniparse(HTMLParser):
	def setID(self, inputID):
		self.myID = inputID
		self.obsolete = 0
		self.deleted = 0
		self.demerged = 0
		self.parsedIDs = []
	def handle_data(self, data):
		if self.demerged:
			matchID = re.search(r'([A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]|[O-Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9])', data.upper())
			if matchID:
				parsedID = matchID.group(0)
				if parsedID.upper() != self.myID:
					self.parsedIDs.append(parsedID)
			return 0
		elif self.obsolete:
			deleted = re.search('deleted',data.lower())
			demerged = re.search('demerged',data.lower())
			if deleted:
				self.deleted = 1
			elif demerged:
				self.demerged = 1
			return 0
		else:
			match = re.search('obsolete',data.lower())
			if match:
				self.obsolete = 1
			return 0
	def returnIDs():
		return parsedIDs

if __name__ == '__main__':
	sys.exit(1)