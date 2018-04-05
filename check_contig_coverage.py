### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###


__usage__ = """
	python check_contig_coverage.py
	--bam <FULL_PATH_TO_INPUT_BAM_FILE>
	--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
	
	optional:
	--bam_is_sorted <PREVENTS_EXTRA_SORTING_OF_BAM_FILE>
	
	bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import sys, os
import matplotlib.pyplot as plt
import numpy as np

# --- end of imports --- #

def analyze_data_and_construct_figure( cov_file, fig_file ):
	"""! @brief load data from coverage file and construct figure """
	
	# --- load read coverage depth of all contigs --- #
	cov_per_contig = {}
	with open( cov_file, "r" ) as f:
		line = f.readline()
		prev_seq = line.split('\t')[0]
		coverage = []
		while line:
			parts = line.strip().split('\t')
			if parts[0] != prev_seq:
				cov_per_contig.update( { prev_seq: np.mean( coverage ) } )
				coverage = []
				prev_seq = parts[0]
			coverage.append( int( parts[2] ) )
			line = f.readline()
		cov_per_contig.update( { prev_seq: np.mean( coverage ) } )
	
	# ---- construct output files --- #
	output_file = fig_file + ".txt"
	with open( output_file, "w" ) as out:
		out.write( "ContigID\tAverageCoverage\n" )
		for contig in sorted( cov_per_contig.keys() ):
			out.write( contig + '\t' + str( cov_per_contig[ contig ] ) + '\n' )
	
	# --- construct figure --- #
	coverage_values = []
	for each in cov_per_contig.values():
		if each <= 1000:
			coverage_values.append( each )
		else:
			print each
	
	fig, ax = plt.subplots()
	
	ax.hist( coverage_values, bins=1000 )
	ax.set_ylabel( "number of contigs" )
	ax.set_xlabel( "average read coverage depth" )
	
	fig.savefig( fig_file, dpi=600 )


def construct_coverage_file( bam_file, output_file, samtools, bedtools, bam_sort_status ):
	
	bam_file = arguments[ arguments.index( '--in' )+1 ]
	output_file = arguments[ arguments.index( '--out' )+1 ]
	
	if '--bam_is_sorted' in arguments:
		sorted_bam_file = bam_file
	else:
		print "sorting BAM file ..."
		sorted_bam_file = output_file + "_sorted.bam"
		cmd = samtools + " sort --threads 8 " + bam_file + " > " + sorted_bam_file
		os.popen( cmd )
	
	# --- calculate read coverage depth per position --- #
	print "calculating coverage per position ...."
	cmd = bedtools + " -d -ibam " + sorted_bam_file + " > " + output_file
	os.popen( cmd )


def main( arguments ):
	"""! @brief run all parts of this script """
	
	bam_file = arguments[ arguments.index('--bam')+1 ]
	output_dir = arguments[ arguments.index('--out')+1 ]
	
	if '--bam_is_sorted' in arguments:
		bam_sort_status = True
	else:
		bam_sort_status = False
	
	samtools = "samtools"	#might need adjustment
	bedtools = "genomeCoverageBed"	#might need adjustment
	
	if output_dir[-1] != '/':
		output_dir += "/"
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	ID = bam_file.split('/')[-1].split('.')[0]
	fig_file = output_dir + ID + "_cov.png"
	cov_file = output_dir + ID + ".cov"
	if not os.path.isfile( cov_file ):
		construct_coverage_file( bam_file, cov_file, samtools, bedtools, bam_sort_status )
	
	analyze_data_and_construct_figure( cov_file, fig_file )


if __name__ == '__main__':
	
	if '--bam' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
