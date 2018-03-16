### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###


__usage__ = """
	python grep_seq_from_fastq.py\n
	--in <FULL_PATH_TO_FASTQ_FILE>
	--out <FULL_PATH_TO_OUTPUT_FASTQ>
	--seq <SEQUENCE_TO_FIND_IN_READ>
	
	bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys

# --- end of imports --- #

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


def main( arguments ):
	"""! @brief runs everything """
	
	fastq_file_in = arguments[ arguments.index('--in')+1 ]
	fastq_file_out = arguments[ arguments.index('--out')+1 ]
	seq = arguments[ arguments.index('--seq')+1 ]
	
	seqs_of_interest = [ seq, revcomp( seq ) ]
	
	with open( fastq_file_out, "w" ) as out:
		with open( fastq_file_in, "r" ) as f:
			line = f.readline()
			while line:
				read_seq = f.readline()
				useless = f.readline()
				qual = f.readline()
				
				status = False
				for seq in seqs_of_interest:
					if seq in read_seq:
						status = True
				if status:
					out.write( line )
					out.write( read_seq )
					out.write( useless )
					out.write( qual )
				
				line = f.readline()


if __name__ == '__main__':
	
	if '--in' in sys.argv and '--out' in sys.argv and '--seq' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
