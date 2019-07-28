### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python blast_filter.py
					--in <FULL_PATH_TO_INPUT_FILE>
					--out <FULL_PATH_TO_OUTPUT_FILE>
					--score <INT, SCORE_CUTOFF>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import sys

# --- end of imports --- #

def main( arguments ):
	"""! @brief filter BLAST(n) results """
	
	input_file = arguments[ arguments.index( '--in' )+1 ]
	output_file = arguments[ arguments.index( '--out' )+1 ]
	score_cutoff = int( arguments[ arguments.index( '--score' )+1 ] )

	with open( output_file, "w" ) as out:
		with open( input_file, "r" ) as f:
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				if float( parts[-1] ) >= score_cutoff:
					out.write( line )
				line = f.readline()


if '--in' in sys.argv and '--out' in sys.argv and '--score' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
