### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
	python get_reads_from_bam.py\n
	--bam <FULL_PATH_TO_BAM>
	--out <FULL_PATH_TO_OUTPUT_DIR>
	optional:
	--min_len <MIN_READ_LENGTH>
	
	bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys

# --- end of imports --- #

def prepare_bam( input_bam, sorted_bam, samtools_path, cpu ):
	"""! @brief extract all reads mapped in pairs and  """
	
	cmd = " ".join( [ samtools_path, "view -f 0x02 -b -u", input_bam, "|", samtools_path, "sort -n --threads", str( cpu ), "-o", sorted_bam, "-" ] )
	print cmd
	os.popen( cmd )


def bam2fastq( sorted_bam, fq1, fq2, bedtools_path ):
	"""! @brief uses bedtools to extract reads from given BAM file """
	
	cmd = " ".join( [ bedtools_path, "bamtofastq -i", sorted_bam, "-fq", fq1, "-fq2", fq2 ] )
	print cmd
	os.popen( cmd )


def select_long_reads( fq1, fq2, fq1_ml, fq2_ml, min_len ):
	"""! @brief get all read pairs which have sufficient length """
	
	with open( fq1_ml, "w", 0 ) as out1:
		with open( fq2_ml, "w", 0 ) as out2:
			with open( fq1, "r" ) as f1:
				with open( fq2, "r" ) as f2:
					line = f1.readline()
					while line:
						s1 = f1.readline()
						u1 = f1.readline()
						q1 = f1.readline()
						h2 = f2.readline()
						s2 = f2.readline()
						u2 = f2.radline()
						q2 = f2.readline()
						if len( s1.strip() ) >= min_len and len( s2.strip() ) >= min_len:
							out1.write( line )
							out1.write( s1 )
							out1.write( u1 )
							out1.write( q1 )
							out2.write( h2 )
							out2.write( s2 )
							out2.write( u2 )
							out2.write( q2 )
						line = f1.readline()


def main( arguments ):
	"""! @brief runs everything """
	
	input_bam = arguments[ arguments.index('--bam')+1 ]
	prefix = arguments[ arguments.index('--out')+1 ]
	
	length_cutoff_status = False
	if '--min_len' in arguments:
		length_cutoff_status = True
		min_len = arguments[ arguments.index('--min_len')+1 ]
	
	cpu = 8
	samtools_path = "samtools"
	bedtools_path = "bedtools"
	
	if prefix[-1] != '/':
		prefix += "/"
	
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	ID = ".".join( input_bam.split('/')[-1].split('.')[:-1] )
	sorted_bam = prefix + ID + ".name_sorted.bam"
	
	prepare_bam( input_bam, sorted_bam, samtools_path, cpu )
	
	fq1 = prefix + ID + "_1.fastq"
	fq2 = prefix + ID + "_2.fastq"
	bam2fastq( sorted_bam, fq1, fq2, bedtools_path )
	
	if length_cutoff_status:
		fq1_ml = prefix + ID + "_long_1.fastq"
		fq2_ml = prefix + ID + "_long_2.fastq"
		select_long_reads( fq1, fq2, fq1_ml, fq2_ml, min_len )



if __name__ == '__main__':
	
	if '--bam' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
