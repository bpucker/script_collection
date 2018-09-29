### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.3 ###

__usage__ = """
						python genome_dot_plot_heatmap.py\n
						--in1 <FASTA1=QUERY>
						--in2 <FASTA2=SUBJECT>
						--out <FULL_PATH_TO_OUTPUT_DIRECTORY>[.]
						--bs <BLOCK_SIZE>
						--qname <NAME_OF_QUERY_WITHOUT_SPACES>
						--sname <NAME_OF_SUBJECT_WITHOUT_SPACES>
						
						OPTIONAL:
						--show	dot plot heatmap will be displayed as interactive figure
						--cite	will not run the script, but display the reference to cite for it\n
						
						bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
						"""

__reference__ = """Pucker 2018"""

import sys, os, re
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from datetime import datetime
from operator import itemgetter

# --- end of imports --- #


def chunk_generator( seq, chunk_size):
    """! @brief generates successive n-sized chunks from input sequence """
    
    for i in xrange( 0, len( seq ), chunk_size):
        yield seq[ i:i + chunk_size]


def construct_seq_block_file( input_file, output_file, block_size ):
	"""! @brief construct file with chunked sequence """
	
	counter = 0
	with open( output_file, "w" ) as out:
		with open( input_file, "r" ) as f:
			f.readline() #remove header
			seq = []
			line = f.readline()
			while line:
				if line[0] == ">":
					seq = "".join( seq )
					chunks = chunk_generator( seq, block_size)
					for each in chunks:
						out.write( '>' + str( counter ).zfill( 12 ) + '\n' + each + '\n' )
						counter += 1
					seq = []
				else:
					seq.append( line.strip() )
				line = f.readline()
			seq = "".join( seq )
			chunks = chunk_generator( seq, block_size)
			for each in chunks:
				out.write( '>' + str( counter ).zfill( 12 ) + '\n' + each + '\n' )
				counter += 1


def load_blast_results( blast_result_file ):
	"""! @brief load blast results for plot """
	
	blast_results = {}
	
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				hit = blast_results[ parts[0] ]
				if float( parts[-1] ) > hit['score']:
					del blast_results[ parts[0] ]
					blast_results.update( { parts[0]: { 'id': parts[0], 'chr': parts[1], 'pos': ( int( parts[8] )+int( parts[9] ) )*0.5, 'score': float( parts[-1] ) } } )
			except KeyError:
				blast_results.update( { parts[0]: { 'id': parts[0], 'chr': parts[1], 'pos': ( int( parts[8] )+int( parts[9] ) )*0.5, 'score': float( parts[-1] ) } } )
			line = f.readline()
	return blast_results


def calculate_color( value ):
	"""! @brief calculates color for given value """
		
	r = 255-int( 255*value**4 )
	g = 255-int( 255*value**4 )
	b = 255
		
	color = '#%02x%02x%02x' % (r, g, b)
	
	return color


def construct_dot_plot( self_blast_results, other_blast_results, output_figure, show_status, block_size, subject_offsets, query_offsets, q_name, s_name ):
	"""! @brief construct dot plot with score weighted points """
	
	x_values = []
	y_values = []
	colors = []
	
	for key in self_blast_results.keys():
		try:
			other_point = other_blast_results[ key ]
			self_point = self_blast_results[ key ]
			score = other_point[ 'score' ] / self_point[ 'score' ]
			point_color = calculate_color( score )
			y = ( int( self_point['id'] ) * block_size )	#query
			x = other_point['pos'] + subject_offsets[ other_point['chr']  ] 	#reference/subject
			
			colors.append( point_color )
			x_values.append( x )
			y_values.append( y )
			
		except KeyError:
			pass
	
	# --- calculation of figure size --- #
	x_len = max( x_values ) - min( x_values )
	y_len = max( y_values ) - min( y_values )
	
	if x_len <= y_len:
		scale1 = x_len / float( y_len )
		scale2 = 1
	else:
		scale1 = 1
		scale2 = y_len / float( x_len )
	
	# --- construction of plot --- #
	fig, ax = plt.subplots( figsize=( 15*scale1, 15*scale2 ) )
	
	ax.set_xlim( min( x_values ), max( x_values ) )
	ax.set_ylim( min( y_values ), max( y_values ) )
	
	#plotting chromosome borders
	for value in query_offsets.values():
		ax.plot( [ min( x_values ), max( x_values ) ], [ value, value ], color="grey", linewidth=1, linestyle="--", alpha=0.1 )	#query chromosome borders
	for value in subject_offsets.values():
		ax.plot( [ value, value ], [ min( y_values ), max( y_values ) ], color="grey", linewidth=1, linestyle="--", alpha=0.1 )	#subject/reference chromosome borders
	
	#plotting data
	ax.scatter( x_values, y_values, c=colors, s=1, marker=".", linewidths=0 )
		
	#adding scale bar
	ax.plot( [ max( x_values )-1000000, max( x_values ) ], [ min( y_values )+1000000, min( y_values )+1000000 ], color="black" )
	ax.text( max( x_values )-1000000, min( y_values )+1200000, "1Mbp", fontsize=2, ha="left", va="bottom" )
	
	ax.set_xlabel( s_name )
	ax.set_ylabel( q_name )
	
	ax.ticklabel_format( axis='y', style="sci", scilimits=(-2,2) )
	ax.ticklabel_format( axis='x', style="sci", scilimits=(-2,2) )
	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	
	patches = []
	for i in range( 11 ):
		patches.append( mpatches.Patch(color=calculate_color( i/10.0 ), label=str(i/10.0)  ) )
	ax.legend( handles=patches, loc='upper left', framealpha=0.5 )	#'lower right'
	
	if show_status:
		plt.show()
	
	fig.savefig( output_figure, dpi=600 )


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


def get_seq_offset( fasta_file ):
	"""! @brief get offset values for all sequences in given file """
	
	offset_values = {}
	with open( fasta_file, "r" ) as f:
		line = f.readline()
		counter = 0
		while line:
			if line[0] == ">":
				offset_values.update( { line.strip()[1:]: counter } )
			else:
				counter += len( line.strip() )
			line = f.readline()
	return offset_values


def main( arguments ):
	"""! @brief runs all parts of this script """
	
	seq_ref_file1 = arguments[  arguments.index( '--in1' )+1 ]	#query
	seq_ref_file2 = arguments[  arguments.index( '--in2' )+1 ]	#subject
	
	prefix = arguments[  arguments.index( '--out' )+1 ]
	if prefix[-1] != "/":
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	if '--show' in arguments:
		show_status = True
	else:
		show_status = False
	
	if '--bs' in arguments:
		block_size = int( arguments[  arguments.index( '--bs' )+1 ] )
	else:
		block_size = 1000
	
	if '--qname' in arguments:
		q_name = arguments[  arguments.index( '--qname' )+1 ] + " [bp]"
	else:
		q_name = "query [bp]"
	
	if '--sname' in arguments:
		s_name = arguments[  arguments.index( '--sname' )+1 ] + " [bp]"
	else:
		s_name = "query [bp]"
	
	output_figure = prefix + "dot_plot_heatmap.png"
	
	# --- splits sequences into chuncks --- #
	seq_block_file1 = prefix + "seq_blocks.fasta"
	if not os.path.isfile( seq_block_file1 ):	
		construct_seq_block_file( seq_ref_file1, seq_block_file1, block_size )
	
	# --- blast DB construction --- #
	self_blast_db = prefix + "self_blast_db"
	other_blast_db = prefix + "other_blast_db"
	os.popen( "makeblastdb -in " + seq_ref_file1 + " -out " +  self_blast_db + " -dbtype nucl" )
	os.popen( "makeblastdb -in " + seq_ref_file2 + " -out " +  other_blast_db + " -dbtype nucl" )
	
	# --- run blast vs. self and vs. other --- #
	self_blast_result_file = prefix + "self_blast_result_file.txt"
	other_blast_result_file = prefix + "other_blast_result_file.txt"
	
	if not os.path.isfile( self_blast_result_file ):	
		os.popen( "blastn -query " + seq_block_file1 + " -db " + self_blast_db + " -out " +  self_blast_result_file + " -outfmt 6 -evalue 0.01 -num_threads 8" )
	if not os.path.isfile( other_blast_result_file ):	
		os.popen( "blastn -query " + seq_block_file1 + " -db " + other_blast_db + " -out " +  other_blast_result_file + " -outfmt 6 -evalue 0.01 -num_threads 8" )
	
	# --- load blast results --- #
	self_blast_results = load_blast_results( self_blast_result_file )
	other_blast_results = load_blast_results( other_blast_result_file )
	
	# --- get offset values for all subject sequences --- #
	query_offsets = get_seq_offset( seq_ref_file1 )
	subject_offsets = get_seq_offset( seq_ref_file2 )
	
	print "Generating dot plot heatmap ..."
	construct_dot_plot( self_blast_results, other_blast_results, output_figure, show_status, block_size, subject_offsets, query_offsets, q_name, s_name )


if __name__ == '__main__':
	
	if '--cite' in sys.argv:
		sys.exit( __reference__ )
		
	elif '--help' in sys.argv or '-h' in sys.argv:
		sys.exit( __usage__ )
	
	elif '--in1' in sys.argv and '--in2' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	
	else:
		sys.exit( __usage__ )
	
	print "all done!"
