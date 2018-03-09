### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.3 ###


__usage__ = """
			python sort_vcf_by_fasta.py\n
			
			--vcf <INPUT_VCF>\n
			--fasta <INPUT_FASTA_FILE>\n
			--output <OUTPUT_VCF_FILE>\n
			
			bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
			"""

import sys
from operator import itemgetter

# --- end of imports --- #

def load_fasta_headers( input_fasta ):
	"""! @brief load fasta headers as they occur within file """
	
	headers = []
	
	with open( input_fasta, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] == '>':
				headers.append( line.strip()[1:].split(' ')[0] )
			line = f.readline()
	return headers


def load_vcf_file( input_vcf ):
	"""! @brief load VCF file into list of dictionaries """
	
	comment_lines = []
	variants = []
	
	with open( input_vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] == '#':
				comment_lines.append( line )
			else:
				parts = line.strip().split('\t')
				try:
					variants.append( { 'chr': parts[0], 'pos': int( parts[1] ), 'line': line } )
				except:
					print line
			line = f.readline()
	return sorted( variants, key=itemgetter( 'chr', 'pos' ) ), comment_lines


def construct_sorted_vcf_file( output_vcf, variants,  headers, comment_lines ):
	"""! @brief construct sorted VCF file """
	
	with open( output_vcf, "w" ) as out:
		for line in comment_lines:
			out.write( line )
		for header in headers:
			for variant in variants:
				if header == variant['chr']:
					out.write( variant['line'] )


def main( parameters ):
	"""! @brief calls all functions for sorting VCF files """
	
	input_fasta = parameters[ parameters.index( '--fasta' )+1 ]
	input_vcf = parameters[ parameters.index( '--vcf' )+1 ]
	output_vcf = parameters[ parameters.index( '--output' )+1 ]
	
	
	variants, comment_lines = load_vcf_file( input_vcf )
	headers = load_fasta_headers( input_fasta )
	construct_sorted_vcf_file( output_vcf, variants,  headers, comment_lines )
	
	

if __name__ == '__main__':
	
	if '--vcf' in sys.argv and '--fasta' in sys.argv and '--output' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
	
