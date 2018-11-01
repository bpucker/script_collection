### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python CRE_detector.py
					--in <FULL_PATH_TO_INPUT_FASTA_FILE>
					--out <FULL_PATH_TO_OUTPUT_TEXT_FILE>
					--cre <FULL_PATH_TO_CRE_FASTA_FILE>
					--mm <NUMBER_OF_PERMITTED_MISSMATCHES [0-2]>
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, re, sys

# --- end of imports --- #


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip().split(" ")[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq } )
					header = line.strip()[1:].split(" ")[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )
	return sequences


def construct_mutated_sequences( seq ):
	"""! @brief construct all possible mutated sequences for a given sequence """
	
	seq = seq.upper()
	
	mutation = [ [ "C", "T", "G" ], [ "A", "T", "G" ], [ "C", "T", "A" ], [ "C", "A", "G" ] ]
	mutation_index = "ACGT"
	
	mutated_sequences = []
	for idx, bp in enumerate( seq ):
		mutation_spectrum = mutation[ mutation_index.index( bp ) ]
		for base in mutation_spectrum:
			try:
				new_seq = seq[:idx] + base + seq[idx+1:]
				mutated_sequences.append( new_seq )
				#here is an option to add the revcomp sequence as well
			except IndexError:
				pass
	#print "number of mutated sequences: " + str( len( mutated_sequences ) ) + " (length of motif: " + str( len( seq ) ) + ")"
	return mutated_sequences


def revcomp( seq ):	#not needed for oriented promotor motifs
	"""! @brief construct reverse complement of sequence """
	
	new_seq = []
	
	bases = { 'a':'t', 't':'a', 'c':'g', 'g':'c', 'n':'n', '-':'-' }
	for nt in seq.lower():
		try:
			new_seq.append( bases[nt] )
		except:
			new_seq.append( 'n' )
	return ''.join( new_seq[::-1] )


def check_promoter_seqs_for_motifs( promoter_seqs, motif_seqs, once_mutated_motifs, twice_mutated_motifs, report_file ):
	"""! @brief check all promoter sequences for all motifs """
	
	promoter_order = promoter_seqs.keys()
	element_order = motif_seqs.keys()
	
	detected_seqs_per_element = []
	for each in element_order:
		detected_seqs_per_element.append( [] )
	
	with open( report_file, "w" ) as out:
		for k, promoter in enumerate( promoter_order ):
			out.write( "analyzing " + promoter + '\n' )
			seq = promoter_seqs[ promoter ]
			for motif in element_order:
				if motif_seqs[ motif ] in seq:
					pos = seq.find( motif_seqs[ motif ] )
					out.write( motif + " detected " + motif + " at position " + str( pos ) + '\n' )
					detected_seqs_per_element[ element_order.index( motif ) ].append( { 'id': promoter.split('_%_%_%_')[0], 'seq': motif_seqs[ motif ], 'pos': pos-len(seq)+50 } )
				else:
					detected_elements = 0
					if len( once_mutated_motifs.keys() ) > 0:
						for each in once_mutated_motifs[ motif ]:
							if each in seq:
								detected_elements_positions = [ m.start() for m in re.finditer( each, seq ) ]
								for pos in detected_elements_positions:
									out.write( motif + " (one mutation) detected " + each + " at position " + str( pos ) + '\n' )
									detected_seqs_per_element[ element_order.index( motif ) ].append( { 'id': promoter.split('_%_%_%_')[0], 'seq': each, 'pos': pos-len(seq)+50 } )
									detected_elements += 1
					if detected_elements:
						if len( twice_mutated_motifs.keys() ) > 0:
							for each in twice_mutated_motifs[ motif ]:
								if each in seq:
									detected_elements_positions = [ m.start() for m in re.finditer( each, seq ) ]
									for pos in detected_elements_positions:
										out.write( motif + " (two mutations) detected " + each + " at position " + str( pos ) + '\n' )
										detected_seqs_per_element[ element_order.index( motif ) ].append( { 'id': promoter.split('_%_%_%_')[0], 'seq': each, 'pos': pos-len(seq)+50 } )
										detected_elements += 1
			out.write( '\n\n\n\n' )
		
		#### THIS MIGHT BE ADDED LATER ###
		# # --- write element output --- #
		# for idx, each in enumerate( element_order ):
			# out.write( each + '\t' + motif_seqs[ each ] + '\n' )
			# for item in detected_seqs_per_element[ idx ]:
				# out.write( "\t".join( map( str, [ item['id'], item['pos'], item['seq'] ] ) ) + '\n' )
			# out.write( '\n\n\n' )


def main( arguments ):
	"""! @brief run all functions """
	
	motif_file = arguments[ arguments.index( '--cre' )+1 ]
	promoter_seq_file = arguments[ arguments.index( '--in' )+1 ]
	report_file = arguments[ arguments.index( '--out' )+1 ]
	if '--mm' in arguments:
		permitted_missmatches = int( arguments[ arguments.index( '--mm' )+1 ] )
	else:
		permitted_missmatches = 1
	
	# --- load all motif sequences --- #
	motif_seqs = load_sequences( motif_file )
	
	# --- mutate all motif sequences once --- #
	once_mutated_versions_per_motif = {}
	if permitted_missmatches > 0:
		for motif in motif_seqs.keys():
			once_mutated_versions_per_motif.update( { motif: construct_mutated_sequences( motif_seqs[ motif ] ) } )
	
	# --- mutate all once mutated motif sequences a second time --- #
	twice_mutated_versions_per_motif = {}
	if permitted_missmatches == 2:
		for motif in motif_seqs.keys():
			mut_seqs = []
			for seq in once_mutated_versions_per_motif[ motif ]:
				mut_seqs += construct_mutated_sequences( seq )
			twice_mutated_versions_per_motif.update( { motif: mut_seqs } )	#contains some revertants as well
	
	# --- load promoter sequences --- #
	promoter_seqs = load_sequences( promoter_seq_file )
	
	# --- check promoter sequences for motifs --- #
	check_promoter_seqs_for_motifs( promoter_seqs, motif_seqs, once_mutated_versions_per_motif, twice_mutated_versions_per_motif, report_file )


if __name__ == '__main__':
	if '--in' in sys.argv and '--out' in sys.argv and '--cre' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
