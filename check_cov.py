### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###

__usage__ = """
					python check_cov.py
					--chr <CHROMOSOME_NAME>
					--start <START_POSITION>
					--end <END_POSITION>
					--cov1 <FULL_PATH_TO_COV_FILE_SPLITTED>
					--cov2 <FULL_PATH_TO_COV_FILE_UNSPLITTED>
					--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import matplotlib.pyplot as plt
import numpy as np
import os, sys

# --- end of imports --- #


def get_cov_region_of_interest( inputfile, outputfile, chromosome, start, end ):
	"""! @brief extract coverage region of interest """
	
	with open( outputfile, "w" ) as out:
		with open( inputfile, "r" ) as f:
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				if parts[0] == chromosome:
					if start <= int( parts[1] ) <= end:
						out.write( line )
				line = f.readline()


def load_cov( cov_file ):
	"""! @brief load coverage from given file """
	
	cov = []
	with open( cov_file, "r" ) as f:
		line = f.readline()
		while line:
			cov.append( int( line.split('\t')[2] ) )
			line = f.readline()
	return cov


def main( arguments ):
	"""! @brief run everything """
	
	chromosome = arguments[ arguments.index( '--chr' )+1 ]
	start = int( arguments[ arguments.index( '--start' )+1 ] )
	end = int( arguments[ arguments.index( '--end' )+1 ] )
	
	output_dir = arguments[ arguments.index( '--out' )+1 ]
	if output_dir[-1] != "/":
		output_dir += "/"
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )

	RNA_cov_file  = arguments[ arguments.index( '--cov1' )+1 ]
	DNA_cov_file = arguments[ arguments.index( '--cov2' )+1 ]

	RNA_selection_file = output_dir + "splitted_selection.cov"
	DNA_selection_file = output_dir + "unsplitted_selection.cov"

	fig_file = output_dir + "coverage_plot.png"
	
	if not os.path.isfile( RNA_selection_file ):
		get_cov_region_of_interest( RNA_cov_file, RNA_selection_file, chromosome, start, end )
	if not os.path.isfile( DNA_selection_file ):
		get_cov_region_of_interest( DNA_cov_file, DNA_selection_file, chromosome, start, end )


	RNA = load_cov( RNA_selection_file )
	DNA = load_cov( DNA_selection_file )

	fig, ax = plt.subplots( figsize=(10,3) )

	ax.plot( np.arange( 0, len( DNA ),1 ), DNA, linestyle="-.", color="black", label="unsplitted" )
	ax.plot( np.arange( 0, len( RNA ),1 ), RNA, linestyle="-.", color="red", label="splitted" )

	ax.set_xlabel( "position in gene [bp]" )

	ax.set_xlim( 0, len( RNA ) )
	ax.set_ylim( 0, max( RNA+DNA ) )

	ax.legend()

	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)

	plt.subplots_adjust( left=0.05, right=0.999, bottom=0.15, top=0.99 )

	fig.savefig( fig_file, dpi=300 )

if '--chr' in sys.argv and '--start' in sys.argv and '--end' in sys.argv and '--cov1' in sys.argv and '--cov2' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
