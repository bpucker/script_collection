### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python sort_gff3_by_pos.py
					--in <FULL_PATH_TO_GFF3_FILE>
					--out <FULL_PATH_TO_GFF3_FILE>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

from operator import itemgetter
import sys

# --- end of imports --- #

def load_lines( input_gff3 ):
	"""! @brief load lines from given GFF3 file """
	
	sort_value = { 'gene': 1, 'transcript': 2, 'mRNA': 2, 'exon': 3, 'CDS': 4, 'protein': 5  }
	lines = []
	with open( input_gff3, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				lines.append( { 'line': line, 'chr': parts[0], 'pos': int( parts[3] ), 'value': sort_value[ parts[2] ] } )
			line = f.readline()
	return sorted( lines, key=itemgetter( 'chr', 'pos', 'value' ) )


def main( parameters ):
	"""! @brief runs everything """

	input_gff3 = parameters[ parameters.index('--in')+1 ]
	output_gff3 = parameters[ parameters.index('--out')+1 ]

	lines = load_lines( input_gff3 )

	with open( output_gff3, "w" ) as out:
		for line in lines:
			out.write( line['line'] )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
