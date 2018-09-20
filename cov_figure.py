### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
	python cov_figure.py
	--cov <FULL_PATH_TO_COVERAGE_FILE>
	--chr <CHROMOSOME_OF_INTEREST>
	--start <START_POSITION>
	--end <END_POSITION>
	
	bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import sys, os, re
import matplotlib.pyplot as plt

# --- end of imports --- #

def load_data( cov_file ):
	"""! @brief load all information from coverage file """
	
	coverage = {}
	with open( cov_file, "r" ) as f:
		line = f.readline()
		seq = line.split('\t')[0]
		cov = []
		while line:
			parts = line.strip().split('\t')
			if parts[0] != seq:
				coverage.update( { seq: cov } )
				seq = parts[0]
				cov = []
			cov.append( float( parts[2] ) )
			line = f.readline()
		coverage.update( { seq: cov } )
	return coverage


def construct_figure( covs, chromosome, start, end, fig_file ):
	"""! @brief construct coverage figure for region of interest """
	
	fig, ax = plt.subplots()
	
	ax.plot( covs, color="green" )
	ax.set_title( chromosome + " " + str( start ) + "-" + str( end ) )
	
	fig.savefig( fig_file, dpi=300 )
	

def main( arguments ):
	"""! @brief run all parts of script """
	
	cov_file = arguments[ arguments.index( '--cov' )+1 ]	#"/vol/agrcourse/members/bpucker/nd1_vs_col.cov"
	
	chromosome = arguments[ arguments.index( '--chr' )+1 ]
	start = int( arguments[ arguments.index( '--start' )+1 ] )
	end = int( arguments[ arguments.index( '--end' )+1 ] )
	
	fig_file = cov_file + chromosome + "_" + str( start ) + "_" + str( end ) + ".png"
	
	coverage = load_data( cov_file )
	construct_figure( coverage[ chromosome ][ start:end ], chromosome, start, end, fig_file )


if __name__ == '__main__':
	
	if '--cov' in sys.argv and '--chr' in sys.argv and '--start' in sys.argv and '--end' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
