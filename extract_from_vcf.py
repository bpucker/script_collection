### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###


__usage__ = """
					python extract_from_vcf.py
					--invcf <FULL_PATH_TO_INPUT_VCF_FILE>
					--outvcf <FULL_PATH_TO_OUTPUT_VCF_FILE>
					--chr <SEQUENCE_NAME>
					--start <START_POSITION>
					--end  <END_POSITION>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import sys

# --- end of imports --- #

def main( parameters ):
	"""! @brief runs everything """

	input_vcf_file = parameters[ parameters.index('--invcf')+1 ]
	output_vcf_file = parameters[ parameters.index('--outvcf')+1 ]

	chromosome = parameters[ parameters.index('--chr')+1 ]
	start = int( parameters[ parameters.index('--start')+1 ] )
	end = int( parameters[ parameters.index('--end')+1 ] )
	
	counter = 0
	with open( output_vcf_file, "w" ) as out:
		with open( input_vcf_file, "r" ) as f:
			line = f.readline()
			while line:
				if line[0] != '#':
					parts = line.strip().split('\t')
					if parts[0] == chromosome:
						if start <= int( parts[1] ) <= end:
							out.write( line )
							counter += 1
				line = f.readline()
	print "number of detected variants: " + str( counter )


if '--invcf' in sys.argv and '--outvcf' in sys.argv and '--chr' in sys.argv and '--start' in sys.argv and '--end' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
