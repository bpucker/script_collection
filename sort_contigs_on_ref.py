### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python sort_contigs_on_ref.py\n
					--contig_file <FULL_PATH_TO_FILE>
					--ref_file <FULL_PATH_TO_FILE>
					--output_dir <FULL_PATH_TO_DIR>
					--species <SOME_SUFFIX_FOR_PSEUDOCHROMOSOME_NAMES>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
"""


import os, sys
from operator import itemgetter

# --- end of imports --- #


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip().split( " " )[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq } )
					header = line.strip()[1:].split( " " )[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )	
	return sequences


def revcomp( seq ):
	"""! @brief construct reverse complement of sequence """
	
	new_seq = []
	
	bases = { 'a':'t', 't':'a', 'c':'g', 'g':'c' }
	for nt in seq.lower():
		try:
			new_seq.append( bases[nt] )
		except:
			new_seq.append( 'n' )
	return ''.join( new_seq[::-1] ).upper()


def blastn_and_sorting( contigs_file, ref_seq_file, prefix, contig_order_file, contigs ):
	""""! @brief run BLASTn and identify optimal order and orientation of contigs """
	
	cmd = "makeblastdb -in " + ref_seq_file + " -out " + prefix+"ref_blastn_db -dbtype 'nucl'"
	os.popen( cmd )
	
	blastn_output_file = prefix + "blastn_output_file.txt"
	cmd = "blastn -query " + contigs_file + " -db " + prefix+"ref_blastn_db -out " + blastn_output_file + " -outfmt 6 -evalue 0.00000001 -num_threads 8" 
	os.popen( cmd )
	
	best_blastn_hits = {}
	
	with open( blastn_output_file, "r" ) as f:
		line = f.readline()
		prev_id = ""
		while line:
			parts = line.strip().split('\t')
			if parts[0] != prev_id:
				if int( parts[8] ) < int( parts[9] ):
					orientation = "fw"
				else:
					orientation = "rv"
				try:
					value = best_blastn_hits[ parts[0] ]['score']
					if value < float( parts[-1] ):
						del best_blastn_hits[ parts[0] ]
						best_blastn_hits.update( { parts[0]: { 'id': parts[0], 'chr': parts[1], 'pos': (int( parts[8] )+int( parts[9] ))/2, 'orientation': orientation, 'score': float( parts[-1] ) } } )
				except KeyError:
					best_blastn_hits.update( { parts[0]: { 'id': parts[0], 'chr': parts[1], 'pos': (int( parts[8] )+int( parts[9] ))/2, 'orientation': orientation, 'score': float( parts[-1] ) } } )
			line = f.readline()
	
	best_blastn_hits = best_blastn_hits.values()	
	sorted_hits = sorted( best_blastn_hits, key=itemgetter( 'chr', 'pos' ) )
	
	with open( contig_order_file, "w" ) as out:
		for hit in sorted_hits:
			try:
				out.write( '\t'.join( map( str, [ hit['id'], hit['chr'], hit['pos'], hit['orientation'], len( contigs[ hit['id'] ] ) ] ) ) + '\n' )
			except:
				print hit['id']


def load_contig_order_file( contig_order_file, assembly_file, pseudo_chromosome_file, species, gap_length=50 ):
	"""! @brief load contig order information from file """
	
	contigs = load_sequences( assembly_file )
	agp_information = []
	
	with open( pseudo_chromosome_file, "w" ) as out:
		with open( contig_order_file, "r" ) as f:
			line = f.readline()
			prev_chr_name = line.split('\t')[1]
			contig_lengths = []
			contig_names = []
			seqs = []
			while line:
				parts = line.strip().split('\t')
				if parts[1] != prev_chr_name:
					agp_information.append( { 'ID': prev_chr_name, 'contig_lengths': contig_lengths, 'contig_names': contig_names } )
					out.write( '>' + prev_chr_name + '_' + species + '\n' + (gap_length*"N").join( seqs ) + '\n' )
					prev_chr_name = parts[1]
					contig_lengths = []
					contig_names = []
					seqs = []
				contig_names.append( parts[0] )
				contig_lengths.append( len( contigs[ parts[0] ] ) )
				if parts[ 3 ] == "fw":
					seqs.append( contigs[ parts[0] ]  )
					
				elif parts[3] == "rv":
					seqs.append( revcomp( contigs[ parts[0] ] ) )
					
				else:
					print "ERROR: SEQ ORIENTATION UNKNOWN!"
				line = f.readline()
			agp_information.append( { 'ID': prev_chr_name, 'contig_lengths': contig_lengths, 'contig_names': contig_names } )
			out.write( '>' + prev_chr_name + '_' + species + '\n' + (gap_length*"N").join( seqs ) + '\n' )
	return agp_information


def write_agp_information_into_file( agp_information, agp_file, gap_length=50 ):
	"""! @brief write the collected scaffolding information into output file """
	
	with open( agp_file, "w" ) as out:
		
		for scaffold in agp_information:
			run_num = 1	#number of element in scaffold
			cur_pos = 1	#start pos of element in scaffold
			con_lengths = scaffold['contig_lengths']
			con_names = scaffold['contig_names']
			
			for idx, con_nam in enumerate( con_names ):
				if idx == ( len( con_names ) - 1 ):	#last element reached
					# --- contig line --- #
					out.write( '\t'.join( [ scaffold['ID'],
								str( cur_pos ),
								str( cur_pos + con_lengths[idx]-1 ),
								str( run_num ),
								"W",
								con_nam,
								"1",
								str( con_lengths[idx] ),
								"+"
								 ]  ) + "\n" )
					
				else:
					# --- contig line --- #
					out.write( '\t'.join( [ scaffold['ID'],
								str( cur_pos ),
								str( cur_pos + con_lengths[idx]-1 ),
								str( run_num ),
								"W",
								con_nam,
								"1",
								str( con_lengths[idx] ),
								"+"
								 ]  ) + "\n" )
					cur_pos += con_lengths[idx]
					run_num += 1
					
					# --- N stretch line --- #
					out.write( '\t'.join( [ scaffold['ID'],
								str( cur_pos ),
								str( cur_pos + gap_length-1 ),
								str( run_num ),
								"N",
								str( gap_length ),
								"scaffold",
								"yes",
								"paired-ends"
								 ]  ) + "\n" )
					cur_pos += gap_length
					run_num += 1


def main( parameters ):
	"""! @brief run all parts of reference-based pseudochromosome construction """
	
	contig_file = parameters[ parameters.index( '--contig_file' )+1 ]
	ref_seq_file = parameters[ parameters.index( '--ref_file' )+1 ]
	
	prefix = parameters[ parameters.index( '--output_dir' )+1 ]
	
	species = parameters[ parameters.index( '--species' )+1 ]
	
	contig_order_file = prefix + "contig_order_TEST.txt"
	
	# --- get positions for all contigs based on reference sequence --- #
	contigs = load_sequences( contig_file )
	blastn_and_sorting( contig_file, ref_seq_file, prefix, contig_order_file, contigs )	
	
	
	# --- construction of AGP file and corresponding multiple fasta file --- #
	pseudo_chromosome_file = prefix + "pseudochromosomes.fasta"
	agp_file = prefix + "final_assembly.agp"
	
	agp_information = load_contig_order_file( contig_order_file, contig_file, pseudo_chromosome_file, species )
	write_agp_information_into_file( agp_information, agp_file )
	

if __name__ == '__main__':
	
	if '--contig_file' in sys.argv and '--ref_file' in sys.argv and '--output_dir' in sys.argv and '--species' in sys.argv:
		main( sys.argv )	
	else:
		sys.exit( __usage__ )
