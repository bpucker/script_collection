### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python get_ZCRs.py
					--ref_cov <FULL_PATH_TO_REFERENCE_COVERAGE_FILE>
					--alt_cov <FULL_PATH_TO_SAMPLE_COVERAGE_FILE>
					--out <FULL_PATH_TO_OUTPUT_FILE>
					"""

import sys

# --- end of imports --- #

def load_coverage( cov_file ):
	"""! @brief load coverage from given file """
	
	coverage = {}
	with open( cov_file, "r" ) as f:
		line = f.readline()
		prev_chr = line.split('\t')[0]
		cov = []
		while line:
			parts = line.strip().split('\t')
			if parts[0] != prev_chr:
				coverage.update( { prev_chr: cov } )
				prev_chr = parts[0]
				cov = []
			cov.append( float( parts[2] ) )
			line = f.readline()
			
		coverage.update( { prev_chr: cov } )
	return coverage


def find_ZCRs( ref_cov, alt_cov, window_size=100 ):
	"""! @brief find zero coverage regions """
	
	ZCRs = []
	for key in alt_cov.keys():
		a_cov = alt_cov[ key ]
		r_cov = ref_cov[ key ]
		
		a_avg = sum( a_cov ) / float( len( a_cov ) )
		r_avg = sum( r_cov ) / float( len( r_cov ) )
		
		print a_avg
		print r_avg
		
		start = 0
		end = start + window_size
		status = True
		
		while status:
			if end > len( r_cov ):
				status = False
				end = len( r_cov )
				if start == end:
					status = False
			
			if sum( a_cov[ start:end ] ) / window_size < 0.05*a_avg:
				in_status = True
				zcr_start = 0 + start
				zcr_end = 0 + start
				while in_status:
					start += window_size
					end += window_size
					zcr_end += window_size
					
					if end > len( r_cov ):
						end = len( r_cov )
						zcr_end = len( r_cov )
						in_status = False
					
					if sum( a_cov[ start:end ] ) / window_size < 0.05*a_avg:
						in_status = True
						if end >= len( r_cov ):
							in_status = False
					else:
						in_status = False
				
				if sum( r_cov[ zcr_start:zcr_end ] ) / ( zcr_end - zcr_start ) > 0.75*r_avg:
					ZCRs.append( { 'chr': key, 'start': zcr_start, 'end': zcr_end, 'ref': sum( r_cov[ zcr_start:zcr_end ] ) / ( zcr_end - zcr_start ), 'alt': sum( a_cov[ zcr_start:zcr_end ] ) / ( zcr_end - zcr_start ) } )
				else:
					print start
				
			start += window_size
			end += window_size
	return ZCRs


def main( arguments ):
	"""! @brief run everything """
	
	ref_cov_file = arguments[ arguments.index( '--ref_cov' )+1 ]
	alt_cov_file = arguments[ arguments.index( '--alt_cov' )+1 ]
	output_file = arguments[ arguments.index( '--out' )+1 ]
	
	ref_cov = load_coverage( ref_cov_file )
	alt_cov = load_coverage( alt_cov_file )


	ZCRs = find_ZCRs( ref_cov, alt_cov )
	print "number of detected ZCRs: " + str( len( ZCRs ) )

	
	with open( output_file, "w" ) as out:
		out.write("Chr\tStart\tEnd\tRefCov\tAltCov\n")
		for ZCR in ZCRs:
			out.write( "\t".join( map( str, [ ZCR['chr'], ZCR['start'], ZCR['end'], ZCR['ref'], ZCR['alt'] ] ) ) + '\n' )


if '--ref_cov' in sys.argv and '--alt_cov' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
