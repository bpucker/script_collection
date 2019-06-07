### Boas Pucker ###
### bpucker@cebitec.uni-bielefelde ###
### v0.1 ###

__usage__ = """
					python transeq.py
					--in <FULL_PATH_TO_INPUT_FILE>
					--out <FULL_PATH_TO_OUTPUT_FILE>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""


import sys

# --- end of imports --- #


def load_genetic_code( genetic_code_file ):
	"""! @brief load genetic code from file """
	
	genetic_code = {}
	
	with open( genetic_code_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			codons = parts[1].split(', ')
			for codon in codons:
				genetic_code.update( { codon: parts[0] } )
			line = f.readline()
	return genetic_code


def load_multiple_fasta_file( fasta_file ):
	"""! @brief load all sequences from multiple fasta file """
	
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


def get_longest_peptide_per_seq( genetic_code, sequences, output_file ):
	"""! @brief get longest encoded peptide per seq and write it into output file """
	
	with open( output_file, "w" ) as out:
		keys = sequences.keys()
		for key in keys:
			seq = sequences[ key ]
			peptide = translate( seq, genetic_code )
			out.write( '>' + key + '\n' + peptide + '\n' )


def translate( seq, genetic_code ):
	"""! @brief translates the given nucleotide sequence into peptide and splits at each star (stop codon) """
	
	seq = seq.upper()
	
	peptide = []
	
	for i in range( int( len( seq ) / 3.0 ) ):
		codon = seq[i*3:i*3+3]
		try:
			peptide.append( genetic_code[ codon ] )
		except:
			peptide.append( "*" )
	return "".join( peptide )


def main( arguments ):
	
	input_file = arguments[ arguments.index( '--in' )+1 ]
	output_file = arguments[ arguments.index( '--out' )+1 ]
	
	genetic_code = {'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'AGT': 'S', 'CAG': 'Q', 'CAA': 'Q', 'CCC': 'P', 'TAG': '*', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CCA': 'P', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S', 'TTA': 'L', 'GAC': 'D', 'TCC': 'S', 'GAA': 'E', 'TCA': 'S', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T', 'TTG': 'L', 'CGT': 'R', 'TAA': '*', 'CGC': 'R'}
	
	sequences = load_multiple_fasta_file( input_file )
	get_longest_peptide_per_seq( genetic_code, sequences, output_file )


if __name__ == '__main__':
	
	if '--in' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
	print "all done!"
