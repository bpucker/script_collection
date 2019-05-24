### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python get_peps_from_gff3.py
					--gff <FULL_PATh_TO_GFF3_FILE>
					--fasta <FULL_PATH_TO_FASTA_FILE>
					--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
					"""

import os, sys
from operator import itemgetter

# --- end of imports --- #


def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip().split(' ')[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq } )
					header = line.strip()[1:].split(' ')[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )
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
				if parts[2] == 'CDS':
					information.append( { 	'chr': parts[0],
															'start': int( parts[3] ),
															'end': int( parts[4] ),
															'orientation': parts[6],
															'parent': parts[-1].split(';')[1]
														} )
			line = f.readline()
	
	# --- sort data by parent --- #
	sorted_data = {}
	for each in information:
		try:
			collected_parts = sorted_data[ each['parent'] ]
			collected_parts.append( each )
			del sorted_data[ each['parent'] ]
			sorted_data.update( { each['parent']: collected_parts } )
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
			seq = ""
			revcomp_status = False
			if transcript[0]['orientation'] == '-':
				revcomp_status = True
			for part in transcript:
				if revcomp_status:
					seq += revcomp( assembly[ part['chr'] ][ part['start']-1:part['end'] ] )
				else:
					seq += assembly[ part['chr'] ][ part['start']-1:part['end'] ]
			out.write( '>' + transcript[0]['parent'].replace("Parent=", "") + '\n' + seq + '\n' )
			

def revcomp( seq ):
	"""! @brief constructs revcomp """
	
	new_seq = []
	dictionary = { 'a':'t', 't':'a', 'c':'g', 'g':'c', 'n':'n' }
	
	for nt in seq.lower():
		new_seq.append( dictionary[ nt ] )
	return ''.join( new_seq[::-1] ).upper()


def main( arguments ):
	"""! @brief extract transcript and pepetides from FASTA and GFF3 """
	
	gff3_input_file = arguments[ arguments.index('--gff')+1 ]
	fasta_input_file = arguments[ arguments.index('--fasta')+1 ]
	output_dir = arguments[ arguments.index('--out')+1 ]
	
	if output_dir[ -1 ] != "/":
		output_dir += "/"
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	pep_output_file = output_dir + "peptides.fasta"
	CDS_output_file = output_dir + "CDS.fasta"
	
	
	seqs = load_sequences( fasta_input_file )
	transcript_information = load_transcript_information_from_gff3( gff3_input_file )
	
	construct_CDS_file( transcript_information, CDS_output_file, seqs )
	
	os.popen( "python transeq.py --in " + CDS_output_file + " --out " + pep_output_file )	#could be integrated into this script
	


if '--gff' in sys.argv and '--fasta' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
