### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
						python dot_plot_heatmap.py\n
						--in <FULL_PATH_TO_ASSEMBLY_FILE>
						--ref <FULL_PATH_TO_REFERENCE_SEQUENCE>
						--out <FULL_PATH_TO_OUTPUT_DIRECTORY>[.]\n
						
						
						OPTIONAL:
						--cite	will not run the script, but display the reference to cite for it\n
						--help	displays this help message
						
						bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
						"""

__reference__ = """Pucker et al., 2018: xxxxx"""

import sys, glob, re, os, time, datetime, shutil
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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
			seq = ""
			line = f.readline()
			while line:
				if line[0] == ">":
					chunks = chunk_generator( seq, block_size)
					for each in chunks:
						out.write( '>' + str( counter ).zfill( 8 ) + '\n' + each + '\n' )
						counter += 1
					seq = ""
				else:
					seq+= line.strip()
				line = f.readline()
			chunks = chunk_generator( seq, block_size)
			for each in chunks:
				out.write( '>' + str( counter ).zfill( 8 ) + '\n' + each + '\n' )
				counter += 1


def load_blast_results( blast_result_file ):
	"""! @brief load blast results for plot """
	
	# --- load BLAST results --- #
	blast_results = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				hits = blast_results[ parts[0] ]
				del blast_results[ parts[0] ]
				hits.append( { 'id': parts[0], 'chr': parts[1], 'pos': ( int( parts[8] )+int( parts[9] ) )*0.5, 'score': float( parts[-1] ) } )
				blast_results.update( { parts[0]: hits } )
			except KeyError:
				blast_results.update( { parts[0]: [ { 'id': parts[0], 'chr': parts[1], 'pos': ( int( parts[8] )+int( parts[9] ) )*0.5, 'score': float( parts[-1] ) } ] } )
			line = f.readline()
	
	# --- process BLAST results --- #
	final_blast_results = {}
	for key in blast_results.keys():
		hits = blast_results[ key ]
		if len( hits ) == 1:
			final_blast_results.update( { key: hits[0] } )
		else:
			hits = sorted( hits, key=itemgetter('score') )
			if hits[-1] > hits[-2]:
				final_blast_results.update( { key: hits[0] } )
			else:
				current = hits[0]
				new_score = current['score']*0.5
				del current['score']
				current.update( { 'score': new_score } )
				final_blast_results.update( { key: current } )
	
	return final_blast_results


def calculate_color( value ):
	"""! @brief calculates color for given value """
		
	r = 255-int( 255*value**4 )
	g = 255-int( 255*value**4 )
	b = 255
		
	color = '#%02x%02x%02x' % (r, g, b)
	
	return color


def  load_seq_lengths( fasta_ref_file ):
	
	sequences = load_sequences( fasta_ref_file )
	lengths = {}
	for key in sequences.keys():
		lengths.update( { key: len( sequences[ key ] ) } )
	return lengths


def construct_dot_plot( blast_results, output_figure, fasta_ref_file ):
	"""! @brief construct dot plot with score weighted points """
		
	chr_lengths = load_seq_lengths( fasta_ref_file )
	chr_names = sorted( chr_lengths.keys() )
	
	# --- construct plot --- #
	fig, ax = plt.subplots( figsize=( 10, 5 ) )
	y_offset = 1.2*( len( chr_names )-1 )
	
	# --- prepare data for plot --- #
	data_x = [ ]
	data_y = [ ]
	scores = []
	for  key in blast_results.keys():
		entry = blast_results[ key ]
		data_x.append( entry['pos']/1000000.0 )
		data_y.append( y_offset-1.2*chr_names.index( entry['chr'] ) )
		scores.append( entry['score'] )
	
	# --- adding chromosomes --- #
	for idx, each in enumerate( chr_names ):
		ax.plot( [ 0,  chr_lengths[ each ]/1000000.0 ], [ y_offset-1.2*idx, y_offset-1.2*idx ] , color="black", linewidth=.5 )
		ax.text( chr_lengths[ each ]/1000000.0 , y_offset-1.2*idx, each, fontsize=5 )
	
	# --- adding genes of interest --- #
	gois = []
	for gene in gois:
		ax.scatter( gene['pos']/1000000.0, y_offset-1.2*chr_names.index( gene['chr'] ), s=1, color="red" )
		ax.text( gene['pos']/1000000.0, y_offset-1.2*chr_names.index( gene['chr'] ), gene['id'], fontsize=5 )
	
	# --- adding variant information --- #
	ax.scatter( data_x, data_y, c=scores, s=1, cmap="bwr" )
	
	# --- general adjustments --- #
	ax.set_xlabel( "chromosome position [Mbp]" )
	
	ax.spines["top"].set_visible(False)
	ax.spines["left"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.set_frame_on(False)
	ax.axes.get_yaxis().set_visible(False)
	
	ax.set_ylim( -0.5, 1.2*len( chr_names )+0.5 )
	
	plt.subplots_adjust( left=-0.02, right=0.98, top=1.0, bottom=0.08 )
	
	fig.savefig( output_figure, dpi=300 )


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


def submit_jobs_to_cluster( prefix, query_file_names, reference_blastn_db, para_jobs, additional_options ):
	"""! @brief submit BLAST jobs for each file to cluster """
	
	IDs_to_check = []
	batch_ID = str( datetime.datetime.now() )[-3:]
	for idx, file_name in enumerate( query_file_names ):
		ID = "B_" + batch_ID + '_' + str( idx ).zfill(4)
		IDs_to_check.append( ID )
		sh_file = prefix + ID + '.sh'
		out_file = prefix + ID + '.out'
		err_file = prefix + ID + '.err'
		
		cmd = "/vol/biotools/bin/blastn -query " + file_name + " -db " + reference_blastn_db + " -out " +  '.'.join( file_name.split('.')[:-1] ) + ".txt " + additional_options
		
		with open( sh_file, "w" ) as out:
				out.write( "#!/bin/bash\n" + " ".join( [ 	"echo " + '"',
																cmd + '"',
																"| qsub -cwd",
																"-N",
																ID,
																"-l vf=1G",
																"-l arch=lx-amd64",
																"-P fair_share",
																"-o",
																out_file,
																"-e",
																err_file
															] ) + '\n'
							)
		os.popen( "chmod +x " + sh_file )
		os.popen( sh_file )
		time.sleep(1)
		os.remove( sh_file )
		waiting_status = True
		while waiting_status:
			qstat = os.popen( "qstat" )
			content = qstat.read()
			qstat_IDs = re.findall( "B_" + batch_ID + "_\d{4}", content )
			counter = 0
			for ID in qstat_IDs:
				if ID in IDs_to_check:
					counter += 1
			if counter < para_jobs:
				waiting_status = False
			else:
				time.sleep( 1 )
	
	waiting_status = True
	while waiting_status:
		qstat = os.popen( "qstat" )
		content = qstat.read()
		qstat_IDs = re.findall( "B_" + batch_ID + "_\d{4}", content )
		waiting_status = False
		for ID in IDs_to_check:
			if ID in qstat_IDs:
				for each in content.split('\n')[2:-1]:
					if ID in each.split()[2] and not 'd' in each.split()[4]:
						waiting_status = True
		time.sleep( 10 )


def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip().split('.')[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq } )
					header = line.strip()[1:].split('.')[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )
	return sequences


def produce_multiple_query_files( query_file, cutoff ):
	"""! @brief produce multiple query files """
	
	prefix = query_file + str( datetime.datetime.now() )[-5:]  + "/"
	os.makedirs( prefix )
	
	sequences = load_sequences( query_file )
	
	query_file_names = []
	
	len_counter = 0
	name_counter = 1
	query_file = prefix + "0".zfill(4) + ".fasta"
	query_file_names.append( query_file )
	out = open( query_file, "w" )
	for idx, seq_id in enumerate( sorted( sequences.keys() ) ):
		if len_counter >= cutoff:
			len_counter = 0
			out.close()
			query_file = prefix + str( name_counter ).zfill(4) + ".fasta"
			query_file_names.append( query_file )
			out = open( query_file, "w" )
			name_counter += 1
		out.write( '>' + seq_id + '\n' + sequences[ seq_id ] + '\n' )
		len_counter += len( sequences[ seq_id ] )
	out.close()
	return prefix, query_file_names


def final_processing( query_file_names, prefix, final_result_file ):
	"""! @brief blt processing of BLAT results for identification of best hit """
	
	result_file_names = []
	for filename in query_file_names:
		result_file_names.append( '.fasta'.join( filename.split('.fasta')[:-1] ) + '.txt' )
	
	cmd1 = "cat " + " ".join( result_file_names ) + " > " + final_result_file
	os.popen( cmd1 )
	
	shutil.rmtree( prefix )


def run_blastn( query_file, prefix, reference_blastn_db, result_file, cluster_status ):
	"""! @brief check inputs and call functions """
	
	if cluster_status:
		para_jobs = 100
		additional_options = ' -word_size 8 -outfmt 6 -evalue 0.01 -max_target_seqs 2'
		
		cutoff=1000000
		
		blast_prefix, query_file_names = produce_multiple_query_files( query_file, cutoff )	#seq length increads compared to BLAT
		
		submit_jobs_to_cluster( blast_prefix, query_file_names, reference_blastn_db, para_jobs, additional_options )
		
		final_processing( query_file_names, blast_prefix, result_file )
	else:
		os.popen( "blastn -query " + query_file + " -db " + reference_blastn_db + " -out " + result_file + additional_options )


def main( arguments ):
	"""! @brief runs all parts of this script """
	
	seq_ref_file1 = arguments[ arguments.index( '--in' )+1 ]
	seq_ref_file2 = arguments[ arguments.index( '--ref' )+1 ]
	
	prefix = arguments[  arguments.index( '--out' )+1 ]
	
	output_figure = prefix + "DOT_PLOT_HEATMAP.png"
	
	block_size = 150
	active = True
	cluster_status = False	#only change for internal use
	
	blast_result_file = prefix + "blast_result_file.txt"
	
	if active:
		
		seq_block_file1 = prefix + "seq_blocks.fasta"
		construct_seq_block_file( seq_ref_file1, seq_block_file1, block_size )
		
		# --- run blast  --- #
		blast_db = prefix + "blast_db"
		os.popen( "makeblastdb -in " + seq_ref_file2 + " -out " +  blast_db + " -dbtype nucl" )
		
		run_blastn( seq_ref_file1, prefix, blast_db, blast_result_file, cluster_status )
		
	## --- load blast results --- #
	blast_results = load_blast_results( blast_result_file )
		
	print "constructing plot ... (this might take a moment)"
	
	#construct plot
	construct_dot_plot( blast_results, output_figure, seq_ref_file2 )


if __name__ == '__main__':
	
	if '--cite' in sys.argv:
		sys.exit( __reference__ )
		
	if '--help' in sys.argv or '-h' in sys.argv:
		sys.exit( __usage__ )
	
	elif '--in' in sys.argv and '--out' in sys.argv and '--ref' in sys.argv:
		main( sys.argv )
	
	else:
		sys.exit( __usage__ )
