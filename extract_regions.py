### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python extract_regions.py
					--in <FULL_PATH_TO_INPUT_FILE>
					--out <FULL_PATH_TO_OUTPUT_FILE>
					--chr <CHROMOSOME_NAME_REGION_OF_INTEREST>
					--start <START_REGION_OF_INTEREST>
					--end <END_REGION_OF_INTEREST>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""


from operator import itemgetter
import sys

# --- end of imports --- #

def load_data( filename ):
	"""! @brief load all data from given file """
	
	data = []
	with open( filename, "r" ) as f:
		header = f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			data.append( { 'chr': parts[0], 'start': int( parts[1] ), 'end': int( parts[2] ), 'rcov': float( parts[3] ), 'acov': float( parts[4] ) } )
			line = f.readline()
	data = sorted( data, key=itemgetter( 'chr', 'start' ) )
	return data, header


def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index( '--in' )+1 ]
	ouput_file = arguments[ arguments.index( '--out' )+1 ]
	chromosome = arguments[ arguments.index( '--chr' )+1 ]
	start = int( arguments[ arguments.index( '--start' )+1 ] )
	end = int( arguments[ arguments.index( '--end' )+1 ] )
		
	data, header = load_data( input_file )
	
	with open( ouput_file, "w" ) as out:
		out.write( header )
		for each in data:
			if each['chr'] == chromosome:
				if each['start'] < end:
					if each['end'] > start:
						out.write( "\t".join( map( str, [ each['chr'], each['start'], each['end'], each['rcov'], each['acov'] ] ) ) + '\n' )


if '--in' in sys.argv and '--out' in sys.argv and '--chr' in sys.argv and '--start' in sys.argv and '--end' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
