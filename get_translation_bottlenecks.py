### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
	python get_translational_bottlenecks.py
	--in <FULL_PATH_TO_INPUT_FILE>
	--codon <FULL_PATH_TO_CODON_USAGE_TABLE>
	--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
						"""

import sys, os
import matplotlib.pyplot as plt

# --- end of imports --- #

def load_all_seqs_from_multiple_fasta_file( filename ):
	"""! @brief load all sequences from multiple fasta file """
	
	data = {}
	
	with open( filename, "r" ) as f:
	 	header = f.readline().strip()[1:].split(' ')[0]
		line = f.readline()
		seq = []
		while line:
			if line[0] == '>':
				data.update( { header: "".join( seq ) } )
				header = line.strip()[1:].split(' ')[0]
				seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		data.update( { header: "".join( seq ) } )
	return data


def analyze_codon_usage( sequences, codon_usage ):
	"""! @brief analyze codon usage in representative CDS """
	
	codon_usage_per_seq = {}
	for key in sequences.keys():
		seq = sequences[ key ].upper()
		if len( seq ) % 3 == 0:
			if not 'N' in seq:
				codon_freq = []
				blocks = [ seq[ i:i+3 ] for i in range( 0, len( seq ), 3 ) ]
				for block in blocks:
					codon_freq.append( codon_usage[ block ] )
				codon_usage_per_seq.update( { key: codon_freq } )
			else:
				print "ERROR: N in coding sequence!"
		else:
			print "ERROR: CDS length is not multiple of 3!"
	return codon_usage_per_seq


def load_codon_usage( codon_usage_file ):
	"""! @brief load codon usage from file """
	
	codon_usage = {}
	with open( codon_usage_file, "r" ) as f:
		f.readline()	#header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			codon_usage.update( { parts[1]: float( parts[3] ) } )
			line = f.readline()
	return codon_usage


def construct_figure( fig_file, ID, values ):
	"""! @brief construct figure to illustrate distribution of rare codons """
	
	fig, ax = plt.subplots( figsize=( 20,5 ) )
	
	ax.set_title( ID )
	for idx, value in enumerate( values ):
		ax.plot( [ idx, idx ], [ 0, value ], color="green" )
	ax.plot( values )
	
	ax.set_xlim( 0.5, len( values )-0.5 )
	plt.subplots_adjust( left=0.03, right=0.98, top=0.99, bottom=0.18, wspace=0.2 )
	
	fig.savefig( fig_file )


def main( arguments ):
	"""! @brief run all functions """
	
	input_seq = arguments[ arguments.index( '--in' )+1 ]
	codon_usage_file = arguments[ arguments.index( '--codon' )+1 ]
	output_dir = arguments[ arguments.index( '--out' )+1 ]
	if output_dir[-1] != '/':
		output_dir += "/"
	if not os.path.exists(output_dir ):
		os.makedirs( output_dir )
	
	codon_usage = load_codon_usage( codon_usage_file )
	sequences = load_all_seqs_from_multiple_fasta_file( input_seq )
	
	codon_usage_per_seq = analyze_codon_usage( sequences, codon_usage )
	output_file = output_dir + "observed_codon_usage.txt"
	with open( output_file, "w" ) as out:
		for ID in codon_usage_per_seq.keys():
			out.write( "\t".join( map( str, [ ID ] + codon_usage_per_seq[ ID ] ) ) + '\n' )
	
	# --- construct figures --- #
	for ID in codon_usage_per_seq.keys():
		fig_file = output_dir + ID + ".png"
		construct_figure( fig_file, ID, codon_usage_per_seq[ ID ] )


if __name__ == '__main__':
	
	if '--in' in sys.argv and '--codon' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
