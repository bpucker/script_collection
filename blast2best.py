### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###

__usage__ = """
					python blast2best.py
					--in <FULL_PATH_TO_INPUT_FILE>
					--out  <FULL_PATH_TO_OUTPUT_FILE>
					"""

import os, sys

# --- end of imports --- #

def load_best_blast_hit( blast_result_file ):
	"""! @brief load best blast hit per query """
	
	best_hits = {}
	
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				data = best_hits[ parts[0] ]
				if float( parts[-1] ) > data['score']:
					del best_hits[ parts[0] ]
					best_hits.update( { parts[0]: { 'score': float( parts[-1] ), 'line': line } } )
			except:
				best_hits.update( { parts[0]: { 'score': float( parts[-1] ), 'line': line } } )
			line = f.readline()
	return best_hits


def main( arguments ):
	
	input_file = arguments[ arguments.index( '--in' ) + 1 ]
	output_file = arguments[ arguments.index( '--out' ) + 1 ]
	
	hits = load_best_blast_hit( input_file )
	with open( output_file, "w" ) as out:
		for key in sorted( hits.keys() ):
			out.write( hits[ key ][ 'line' ] )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
