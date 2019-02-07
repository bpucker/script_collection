### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
	python tree.py
	--in <FULL_PATH_TO_INPUT_FILE>
	--out <FULL_PATH_TO_OUTPUT_DIR>
	
	optional:
	--occ <FLOAT, occupancy required per alignment column>
	--name <STRING, prefix for final alignment file>
					"""

import os, sys

# --- end of imports --- #


def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip().replace( " ", "" ).replace('(', '[').replace(')', ']')
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq } )
					header = line.strip()[1:].replace( " ", "" ).replace('(', '[').replace(')', ']')
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )
	return sequences


def modify_FASTA_names( input_file, output_file, mapping_table_file ):
	"""! @brief modify names in FASTA file """
	
	seqs = load_sequences( input_file )
	mapping_table = {}
	with open( output_file, "w" ) as out:
		for idx, key in enumerate( seqs.keys() ):
			out.write( '>seq' + str( idx ).zfill(5) + '\n' + seqs[ key ] + '\n' )
			mapping_table.update( { key: 'seq' + str( idx ).zfill(5) } )
	
	with open( mapping_table_file, "w" ) as out:
		for key in mapping_table.keys():
			out.write( key + '\t' + mapping_table[ key ] + '\n' )


def load_mapping_table( mapping_table_file ):
	"""! @brief load mapping table """
	
	mapping_table = {}
	with open( mapping_table_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			mapping_table.update( { parts[1]: parts[0] } )
			line = f.readline()
	return mapping_table


def modify_names_in_tree( input_tree_file, output_tree_file, mapping_table_file ):
	"""! @brief replace short names in tree by original names """
	
	mapping_table = load_mapping_table( mapping_table_file )
	
	with open( input_tree_file, "r" ) as f:
		tree = f.read()
	
	for key in mapping_table.keys():
		tree = tree.replace( key, mapping_table[ key ] )
	
	with open( output_tree_file, "w" ) as out:
		out.write( tree )


def main( arguments ):
	"""! @brief handle everything """
	
	input_file = arguments[ arguments.index('--in') + 1 ]
	output_dir = arguments[ arguments.index('--out') + 1 ]
	
	mafft = "mafft"
	pxclsq = "pxclsq"
	fasttree = "FastTree"
	
	if '--occ' in arguments:
		occupancy = float( arguments[ arguments.index('--occ') + 1 ] )
	else:
		occupancy = 0.3
	
	if '--name' in arguments:
		name = arguments[ arguments.index('--name') + 1 ]
	else:
		name = ""
	
	if output_dir[-1] != '/':
		output_dir += "/"
	
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	mod_FASTA = output_dir + "names_modified.fasta"
	mapping_table_file = output_dir + "seq_names_mapping_table.txt"
	modify_FASTA_names( input_file, mod_FASTA, mapping_table_file )
	
	alignment_file = mod_FASTA + ".aln"
	os.popen( " ".join( [ mafft, mod_FASTA, ">", alignment_file ] ) )
	
	clean_alignment_file = alignment_file + ".cln"
	os.popen( " ".join( [ pxclsq, "-s", alignment_file, "-o", clean_alignment_file, "-p", str( occupancy ) ] ) )
	
	tree_file = clean_alignment_file + ".tre"
	os.popen( " ".join( [ fasttree, "-wag -nosupport <", clean_alignment_file, ">", tree_file ] ) )
	
	output_tree_file = output_dir + name + "FINAL_TREE.tre"
	modify_names_in_tree( tree_file, output_tree_file, mapping_table_file )


if __name__ == '__main__':
	
	if '--in' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
