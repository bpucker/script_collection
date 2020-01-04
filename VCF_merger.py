### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###

__usage__ = """
					python VCF_merger.py
					--in <INPUT_FOLDER>
					--out <OUTPUT_FILE>
					
					"""

import os, glob, sys

# --- end of imports --- #


def load_variants( filename ):
	"""! @brief load all variants from given file """
	
	tmp_variants = {}
	with open( filename, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				tmp_variants.update( { parts[0] + "_%_" + parts[1].zfill( 9 ): line } )
			line = f.readline()
	return tmp_variants


def main( arguments ):
	
	vcf_folder = arguments[ arguments.index( '--in' )+1 ]
	output_vcf = arguments[ arguments.index( '--out' )+1 ]
	
	if vcf_folder[-1] != '/':
		vcf_folder += "/"

	vcfs = glob.glob( vcf_folder + "*.vcf" )

	print "number of VCF files: " + str( len( vcfs ) )

	# --- load variants from all files --- #
	variants = {}
	for filename in vcfs:
		variants.update( load_variants( filename ) )

	with open( output_vcf, "w" ) as out:
		for key in sorted( variants.keys() ):
			out.write( variants[ key ] )

	print "total number of final variants: " + str( len( variants.keys() ) )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
