### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
	python candidate_gene_identification.py\n
	--query <FULL_PATH_TO_QUERY_FILE>
	--pep <FULL_PATH_TO_SUBJECT_PEPTIDE_FILE>
	--prefix <FULL_PATH_TO_OUTPUT_DIRECTORY>
	
	bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import sys, os
from operator import itemgetter

# --- end of imports --- #

def get_all_subject_IDs( blast_result_file ):
	"""! @brief get all subject IDs """
	
	all_subjects = {}
	
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				all_subjects[ parts[1] ]
			except KeyError:
				all_subjects.update( {  parts[1]: None} )
			line = f.readline()
	return all_subjects


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


def construct_pre_alignment_file( output_file, peps, queries, subject_IDs ):
	"""! @brief write all sequences into output file to enable multiple alignment via MAFFT """
	
	with open( output_file, "w" ) as out:
		for key in queries.keys():
			out.write( '>' + key + '\n' + queries[ key ] + '\n' )
		for key in subject_IDs:
			out.write( '>' + key + '\n' + peps[ key ] + '\n' )


def load_blast_results( result_file ):
	"""! @brief load detailed informaiton about BLASTp results """
	
	# --- loading BLAST results --- #
	blast_results = {}
	with open( result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				value = blast_results[ parts[1] ]['score']
				if float( parts[-1] ) > value:
					del blast_results[ parts[1] ]
					blast_results.update( { parts[1]: { 'query': parts[0], 'id': parts[1], 'score': float( parts[-1] ) } } )
			except KeyError:
				blast_results.update( { parts[1]: { 'query': parts[0], 'id': parts[1], 'score': float( parts[-1] ) } } )
			line = f.readline()
	
	# --- get best BLAST results --- #
	final_blast_results = {}
	for key in blast_results.keys():
		try:
			current_value = final_blast_results[ key ]['score']
			if current_value < blast_results[ key ][ 'score' ]:
				del final_blast_results[ key ]
				final_blast_results.update( { key: blast_results[ key ] } )
		except KeyError:
			final_blast_results.update( { key: blast_results[ key ] } )
	return final_blast_results


def main( arguments ):
	"""! @brief run everything """
	
	query_file = arguments[ arguments.index( '--query' )+1 ]
	pep_file = arguments[ arguments.index( '--pep' )+1 ]
	prefix = arguments[ arguments.index( '--prefix' )+1 ]
	
	# --- these paths could be adjusted to enable adaptation to other systems --- #
	makeblast_db_path = "makeblastdb"
	blastp_path = "blastp"
	mafft_path = "mafft"
	pxclsq_path = "pxclsq"
	fast_tree_path = "FastTree"
	
	pre_alignment_file = prefix + "seqs_for_alignment.fasta"
	alignment_file = prefix + "aligned_seqs.fasta"
	clean_alignment_file = prefix + "aligned_seqs.fasta.clean"
	tree_file = prefix + "tree_of_seqs.tree"
	blastp_based_classification_file = prefix + "BLASTp_based_classification.txt"
	
	# ---- constructing BLAStp database --- #
	blastp_db = prefix + "blastp_db"
	os.popen( makeblast_db_path + " -in " + pep_file + " -out " + blastp_db + " -dbtype prot" )
	
	# --- running BLASTp against query file --- #
	blast_result_file = prefix + "blast_results.txt"
	os.popen( blastp_path + " -query " + query_file + " -db " + blastp_db + " -out " + blast_result_file + " -outfmt 6 -evalue 0.00001" )
	
	subject_IDs = get_all_subject_IDs( blast_result_file )
	peps = load_all_seqs_from_multiple_fasta_file( pep_file )
	queries = load_all_seqs_from_multiple_fasta_file( query_file )
	
	# --- identify best hit per query and write it to file --- #
	final_blast_results = load_blast_results( blast_result_file )
	sorted_hits = sorted( final_blast_results.values(), key=itemgetter( 'score', 'query' ) )[ ::-1 ]
	with open( blastp_based_classification_file, "w" ) as out:
		out.write( "SubjectID\tQueryID\tBLAST_Score\n" )
		for hit in sorted_hits:
			out.write( hit['id'] + '\t' + hit['query'] + '\t' + str( hit['score'] ) + '\n'  )
	
	# --- construct phylogenetic tree to confirm BLASTp-based classification --- #
	if len( subject_IDs ) + len( queries.keys() ) > 2:
		os.chdir( prefix )
		# --- prepare file with all sequences to enable aignment --- #
		construct_pre_alignment_file( pre_alignment_file, peps, queries, subject_IDs )
		
		# --- constructing MAFFT alignment --- #
		os.popen( mafft_path + " " + pre_alignment_file + " > " + alignment_file )
		
		# --- removing columns with almost no information --- #
		os.popen( pxclsq_path + " -s " + alignment_file + " -o " + clean_alignment_file + " -p 0.1" )
		
		# --- constructing phylogenetic tree --- #
		os.popen( fast_tree_path + " -wag -nosupport < " + clean_alignment_file + " > " + tree_file )
	else:
		print "ERROR: insufficient number of sequences to construct a tree"


if __name__ == '__main__':
	
	if '--query' in sys.argv and '--pep' in sys.argv and '--prefix' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
