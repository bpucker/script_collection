### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

import sys, gzip

__usage__ = """
		python split_FASTQ.py\n
		--in_file <FULL_PATH_TO_FILE>
		
		splits FASTQ files with interlocked mates into separate files for FW and RV reads
		
		bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

def split_compressed_file( input_file, out_fw_file, out_rv_file ):
	"""! @brief split compressed FASTQ file into files for fw and rv reads, respectively """
	
	with gzip.open( out_fw_file, "wb" ) as fw:
		with gzip.open( out_rv_file, "wb" ) as rv:
			with gzip.open( input_file, "rb" ) as f:
				line = f.readline()
				while line:
					# --- fw read --- #
					fw.write( line )
					fw.write( f.readline() )
					fw.write( f.readline() )
					fw.write( f.readline() )
					
					# --- rv read --- #
					rv.write( f.readline() )
					rv.write( f.readline() )
					rv.write( f.readline() )
					rv.write( f.readline() )
					
					line = f.readline()


def split_uncompressed_file( input_file, out_fw_file, out_rv_file ):
	"""! @brief split compressed FASTQ file into files for fw and rv reads, respectively """
	
	with open( out_fw_file, "w" ) as fw:
		with open( out_rv_file, "w" ) as rv:
			with open( input_file, "r" ) as f:
				line = f.readline()
				while line:
					# --- fw read --- #
					fw.write( line )
					fw.write( f.readline() )
					fw.write( f.readline() )
					fw.write( f.readline() )
					
					# --- rv read --- #
					rv.write( f.readline() )
					rv.write( f.readline() )
					rv.write( f.readline() )
					rv.write( f.readline() )
					
					line = f.readline()


def main( arguments ):
	"""! @brief runs everything """
	
	input_file = arguments[ arguments.index( '--in_file' )+1 ]
	directory = "/".join( input_file.split('/')[:-1] ) + "/"
	ID = input_file.split('/')[-1].split('.')[0]
	if 'gz' in input_file.split('.'):
		compressed = True
	elif 'gzip' in input_file.split('.'):
		compressed = True
	else:
		compressed = False
	
	if compressed:
		out_fw_file = directory + ID + "_1.fastq.gz"
		out_rv_file = directory + ID + "_2.fastq.gz"
		split_compressed_file( input_file, out_fw_file, out_rv_file )
	else:
		out_fw_file = directory + ID + "_1.fastq"
		out_rv_file = directory + ID + "_2.fastq"
		split_uncompressed_file( input_file, out_fw_file, out_rv_file )



if __name__ == '__main__':
	
	if '--in_file' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
