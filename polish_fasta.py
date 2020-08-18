### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python polish_fasta.py
					--in <INPUT_FASTA_FILE>
					--out <OUTPUT_FASTA_FILE>
					"""


import os, sys

# --- end of imports --- #

def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	names = []
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip()
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
				sequences.update( { header: "".join( seq ) } )
				names.append( header )
				header = line.strip()[1:]
				seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )
		names.append( header )
	return sequences, names


def main( arguments ):
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	
	seqs, names = load_sequences( input_file )
	
	with open( output_file, "w" ) as out:
		for name in names:
			out.write( '>' + name + '\n' + seqs[ name ] + "\n" )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
