### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.4 ###

__usage__ = """
					python get_peps_from_gff3.py
					--gff <FULL_PATh_TO_GFF3_FILE>
					--fasta <FULL_PATH_TO_FASTA_FILE>
					--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
					
					optional:
					--name <NAME_OF_OUTPUT>
					
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys
from operator import itemgetter

# --- end of imports --- #


def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip()
		if " " in header:
			header = header.split(' ')[0]
			if "\t" in header:
				header = header.split('\t')[0]
		elif "\t" in header:
				header = header.split('\t')[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ).upper() } )
					header = line.strip()[1:]
					if " " in header:
						header = header.split(' ')[0]
						if "\t" in header:
							header = header.split('\t')[0]
					elif "\t" in header:
						header = header.split('\t')[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ).upper() } )
	return sequences


def load_transcript_information_from_gff3( gff3_input_file ):
	"""! @brief load all transcript information from gff3 file """
	
	# --- load all data from file --- #
	information = []
	with open( gff3_input_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if len( parts ) > 5:
					if parts[2].upper() == 'CDS':
						if len( parts[-1] ) > len( "Parent=" ):
							if ";" in parts[-1]:
								parent = False
								subparts = parts[-1].split(';')
								for subp in subparts:
									if "Parent=" in subp:
										parent = subp.replace( "Parent=", "" )
									elif "transcript_id " in subp:
										parent = subp.split('"')[1]
									
								if parent:
									information.append( { 	'chr': parts[0],
																		'start': int( parts[3] ),
																		'end': int( parts[4] ),
																		'orientation': parts[6],
																		'parent': parent
																	} )
								else:
									print "no parent detected - " + line
							else:
								print "only one field - " + line
							
			line = f.readline()
	
	# --- sort data by parent --- #
	sorted_data = {}
	for each in information:
		try:
			sorted_data[ each['parent'] ].append( each )
		except KeyError:
			sorted_data.update( { each['parent']: [ each ] } )
	
	final_data = []
	for key in sorted_data.keys():
		if sorted_data[ key ][0] ['orientation'] == '+':
			final_data.append( sorted( sorted_data[ key ], key=itemgetter('start') ) )
		else:
			final_data.append( sorted( sorted_data[ key ], key=itemgetter('start') )[::-1] )
	return final_data


def construct_CDS_file( transcript_info, CDS_file, assembly ):
	"""! @brief construct file with all sequences for translation """
	
	with open( CDS_file, "w" ) as out:
		for transcript in transcript_info:
			seq = []
			revcomp_status = False
			if transcript[0]['orientation'] == '-':
				revcomp_status = True
			for part in transcript:
				if revcomp_status:
					seq.append( revcomp( assembly[ part['chr'] ][ part['start']-1:part['end'] ] ) )
				else:
					seq.append( assembly[ part['chr'] ][ part['start']-1:part['end'] ] )
			out.write( '>' + transcript[0]['parent'].replace("Parent=", "") + '\n' + "".join( seq ) + '\n' )


def revcomp( seq ):
	"""! @brief constructs revcomp """
	
	new_seq = []
	dictionary = { 'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N' }
	
	for nt in seq:
		try:
			new_seq.append( dictionary[ nt ] )
		except KeyError:
			print "ERROR: " + nt
			new_seq.append( "N")
	return ''.join( new_seq[::-1] )


def translate( seq, genetic_code ):
	"""! @brief translates the given nucleotide sequence into peptide and splits at each star (stop codon) """
	
	seq = seq.upper()
	peptide = []
	for i in range( int( len( seq ) / 3.0 ) ):
		codon = seq[i*3:i*3+3]
		try:
			peptide.append( genetic_code[ codon ] )
		except:
			peptide.append( "*" )
	peptide =  "".join( peptide )
	if sum( [ peptide[0] != "M", "*" in peptide[:-1] ] ) > 0:
		peptide2 = []
		for i in range( int( ( len( seq )-1 ) / 3.0 ) ):
			codon = seq[1+i*3:1+i*3+3]
			try:
				peptide2.append( genetic_code[ codon ] )
			except:
				peptide2.append( "*" )
		peptide2 =  "".join( peptide2 )
		if sum( [ peptide2[0] != "M", "*" in peptide2[:-1] ] ) > 0:
			peptide3 = []
			for i in range( int( ( len( seq )-2 ) / 3.0 ) ):
				codon = seq[2+i*3:2+i*3+3]
				try:
					peptide3.append( genetic_code[ codon ] )
				except:
					peptide3.append( "*" )
			peptide3 =  "".join( peptide3 )
			
			pep_options = []
			if '*' in peptide:
				pep_options.append( { 'seq': peptide.split('*')[0], 'stopps': peptide.count('*'), 'len': len( peptide.split('*')[0] ) } )
			else:
				pep_options.append( { 'seq': peptide, 'stopps': peptide.count('*'), 'len': len( peptide ) } )
			if "*" in peptide2:
				pep_options.append( { 'seq': peptide2.split('*')[0], 'stopps': peptide2.count('*'), 'len': len( peptide2.split('*')[0] ) } )
			else:
				pep_options.append( { 'seq': peptide2, 'stopps': peptide2.count('*'), 'len': len( peptide2 ) } )
			if "*" in peptide3:
				pep_options.append( { 'seq': peptide3.split('*')[0], 'stopps': peptide3.count('*'), 'len': len( peptide3.split('*')[0] ) } )
			else:
				pep_options.append( { 'seq': peptide3, 'stopps': peptide3.count('*'), 'len': len( peptide3 ) } )
			return sorted( pep_options, key=itemgetter( 'len' ) )[-1]['seq']
		else:
			pep_options = []
			if '*' in peptide:
				pep_options.append( { 'seq': peptide.split('*')[0], 'stopps': peptide.count('*'), 'len': len( peptide.split('*')[0] ) } )
			else:
				pep_options.append( { 'seq': peptide, 'stopps': peptide.count('*'), 'len': len( peptide ) } )
			if "*" in peptide2:
				pep_options.append( { 'seq': peptide2.split('*')[0], 'stopps': peptide2.count('*'), 'len': len( peptide2.split('*')[0] ) } )
			else:
				pep_options.append( { 'seq': peptide2, 'stopps': peptide2.count('*'), 'len': len( peptide2 ) } )
			return sorted( pep_options, key=itemgetter( 'len' ) )[-1]['seq']
	else:
		return peptide


def transeq( input_file, output_file, strict_start, strict_end ):
	"""! @brief run translation of coding sequences """	
	genetic_code = {'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'AGT': 'S', 'CAG': 'Q', 'CAA': 'Q', 'CCC': 'P', 'TAG': '*', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CCA': 'P', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S', 'TTA': 'L', 'GAC': 'D', 'TCC': 'S', 'GAA': 'E', 'TCA': 'S', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T', 'TTG': 'L', 'CGT': 'R', 'TAA': '*', 'CGC': 'R'}
	
	sequences = load_sequences( input_file )
	with open( output_file, "w" ) as out:
		keys = sequences.keys()
		for key in keys:
			seq = sequences[ key ]
			if len( seq ) > 9:
				peptide = translate( seq, genetic_code )
				out.write( '>' + key + '\n' + peptide + '\n' )
			else:
				print key + " - too short!"


def main( arguments ):
	"""! @brief extract transcript and pepetides from FASTA and GFF3 """
	
	gff3_input_file = arguments[ arguments.index('--gff')+1 ]
	fasta_input_file = arguments[ arguments.index('--fasta')+1 ]
	output_dir = arguments[ arguments.index('--out')+1 ]
	
	if '--name' in arguments:
		name = arguments[ arguments.index('--name')+1 ]
	else:
		name = "xxx"
	
	if output_dir[ -1 ] != "/":
		output_dir += "/"
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	pep_output_file = output_dir + name + ".pep.fasta"
	CDS_output_file = output_dir + name + ".cds.fasta"
	
	
	seqs = load_sequences( fasta_input_file )
	transcript_information = load_transcript_information_from_gff3( gff3_input_file )
	
	strict_start = False	#option not implemented yet
	strict_end = False	#option not implemented yet
	
	construct_CDS_file( transcript_information, CDS_output_file, seqs )
	transeq( CDS_output_file, pep_output_file, strict_start, strict_end )


if '--gff' in sys.argv and '--fasta' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
