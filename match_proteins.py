### Boas Pucker ###
### v0.3 ###
### bpucker@cebitec.uni-bielefeld.de ###

__usage__ = """ python match_proteins.py\n
				
				--prefix <FULL_PATH_TO_DIRECTORY_FOR_TMP_DATA_AND_RESULTS>\n
				--input1 <FULL_PATH_TO_INPUT_FILE1>\n
				--input2 <FULL_PATH_TO_INPUT_FILE2>
				
				feature requests and bug reports: bpucker@cebitec.uni-bielefeld.de
			"""

import re, os, sys
from operator import itemgetter


# --- end of imports --- #

def load_results_from_BLAST_result_file( BLAST_result_file, cutoff=0.9999 ):
	"""! @brief load data from BLAST result file """
	
	data = {}
	
	with open( BLAST_result_file, "r" ) as f:
		line = f.readline()
		prev_query = line.split('\t')[0]
		hits = []
		while line:
			parts = line.strip().split('\t')
			if parts[0] != prev_query:
				sorted_hits = sorted( hits, key=itemgetter( 'score' ) )
				if len( sorted_hits ) > 1:
					if ( sorted_hits[-2]['score'] / sorted_hits[-1]['score'] ) < cutoff:
						data.update( { sorted_hits[-1]['query']: sorted_hits[-1]['subject'] } )
				else:
					data.update( { sorted_hits[-1]['query']: sorted_hits[-1]['subject'] } )
				hits = []
				prev_query = parts[0]
			hits.append( { 'query': parts[0], 'subject': parts[1], 'score': float( parts[-1] ) } )
			line = f.readline()
		sorted_hits = sorted( hits, key=itemgetter( 'score' ) )
		if len( sorted_hits ) > 1:
			if ( sorted_hits[-2]['score'] / sorted_hits[-1]['score'] ) > cutoff:
				data.update( { sorted_hits[-1]['query']: sorted_hits[-1]['subject'] } )
		else:
			data.update( { sorted_hits[-1]['query']: sorted_hits[-1]['subject'] } )
	#print "entries in data: " + str( len( data.keys() ) )
	return data


def compare_datasets( data1, data2, outputfile, best_score ):
	"""! @brief compares datasets and identifies bidirectional best hits """
	
	seq_IDs_of_interest = []
	
	counter = 0
	keys = data1.keys()
	rbhs = {}
	with open( outputfile, "w" ) as out:
		out.write( "ID1\tID2\tstatus\tscore\n" )
		# --- identify RBHs --- #
		for key in keys:	#key=candidate gene
			try:
				value = data1[ key ]	#value=contig_ID
				try:
					other_value = data2[ value ]	#other_value=candidate_gene_ID
					if key == other_value:
						counter += 1
						out.write( key + '\t' + value + '\tRBH\t' + str( best_score[ key ] ) + '\n' )
						rbhs.update( { key: None } )
						seq_IDs_of_interest.append( value )
				except:
					pass
			except:
				pass
		#print "number of RBH matches: " + str( counter )
		
		# --- identify additional matches --- #
		for key in keys:
			try:
				rbhs[ key ]
			except KeyError:
				out.write( key + '\t' + data1[ key ] + '\tBBH\t' + str( best_score[ key ] ) + '\n' )
				counter += 1
		#print "final number of all matches: " + str( counter )
	return seq_IDs_of_interest


def load_multiple_fasta_file( fasta_file ):
	"""!@brief load content of multiple fasta file """
	
	content = {}
	with open( fasta_file, "r" ) as f:
		header = f.readline().strip()[1:].split(' ')[0]
		line = f.readline()
		seq = ""
		while line:
			if line[0] == '>':
				content.update( { header: seq } )
				header = line.strip()[1:].split(' ')[0]
				seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		content.update( { header: seq } )
	return content


def  load_best_hit_score( seq_file1_blast_result_file ):
	"""! @brief load best score per hit """
	
	best_score = {}
	with open( seq_file1_blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				if best_score[ parts[0] ] < float( parts[-1] ):
					best_score[ parts[0] ] = float( parts[-1] )
			except KeyError:
				best_score.update( { parts[0]: float( parts[-1] ) } )
			line = f.readline()
	return best_score


def identify_protein_matches( parameters ):
	"""! @brief identifies RBHs between given data sets """
	
	prefix = parameters[ parameters.index( '--prefix' )+1 ]
	if prefix[-1] != '/':
		prefix += '/'
	
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	seq_file1 = parameters[ parameters.index( '--input1' )+1 ]
	seq_file2 = parameters[ parameters.index( '--input2' )+1 ]
	
	if not os.path.isfile( seq_file1 ):
		sys.exit( "ERROR: input file1 not detected!" )
	if not os.path.isfile( seq_file2 ):
		sys.exit( "ERROR: input file2 not detected!" )
	
	if '--cpu' in parameters:
		cpu = int( parameters[ parameters.index( '--cpu' )+1 ] )
	else:
		cpu = 8
	
	RBH_file = prefix + "RBH_file.txt"
	
	seq_file1_db = prefix + "seq_file1_db"
	seq_file2_db= prefix + "seq_file2_db"
	
	seq_file1_blast_result_file = prefix + "seq_file1_blast_result_file.txt"
	seq_file2_blast_result_file = prefix + "seq_file2_blast_result_file.txt"
	
	# --- identify RBHs --- #
	#print "db construction ... please wait!"
	db_command1 = "/vol/biotools/bin/makeblastdb -in " + seq_file1 + " -out " + seq_file1_db + " -dbtype 'prot' -parse_seqids"
	db_command2 = "/vol/biotools/bin/makeblastdb -in " + seq_file2 + " -out " + seq_file2_db + " -dbtype 'prot' -parse_seqids"
	
	os.popen( db_command1 )
	os.popen( db_command2 )
	#print " ... done."
	
	
	#print "BLAST1 ... "
	BLAST_command1 = "/vol/biotools/bin/blastp -query " + seq_file1 + " -db " + seq_file2_db + " -out " + seq_file1_blast_result_file + " -outfmt 6 -evalue 0.0001 -num_threads " + str( cpu )
	os.popen( BLAST_command1 )
	
	#print "BLAST2 ... "
	BLAST_command2 = "/vol/biotools/bin/blastp -query " + seq_file2 + " -db " + seq_file1_db + " -out " + seq_file2_blast_result_file + " -outfmt 6 -evalue 0.0001 -num_threads " + str( cpu )
	os.popen( BLAST_command2 )
	
	
	#print "analyzing BLAST results ... please wait!"
	data1 = load_results_from_BLAST_result_file( seq_file1_blast_result_file )
	data2 = load_results_from_BLAST_result_file( seq_file2_blast_result_file )
	best_score = load_best_hit_score( seq_file1_blast_result_file )
	seq_IDs_of_interest = compare_datasets( data1, data2, RBH_file, best_score )
	

if __name__ == '__main__':
	
	if '--prefix' in sys.argv and '--input1' in sys.argv and '--input2' in sys.argv:
		identify_protein_matches( sys.argv )
	else:
		sys.exit( __usage__ )
	
	#print "all done!"
