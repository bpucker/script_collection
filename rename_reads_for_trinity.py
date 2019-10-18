### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
				python rename_reads_for_trinity.py
				--fwin <FW_READ_INPUT_FILENAME>
				--rvin <RV_READ_INPUT_FILENAME>
				--fwout <FW_READ_OUTPUT_FILENAME>
				--rvout <RV_READ_OUTPUT_FILENAME>
				"""

import sys, gzip

# --- end of imports --- #

def main( arguments ):
	"""! @brief generate read names for fw and rv read based on given read number in file """
	
	input_fw_file = arguments[ arguments.index( '--fwin' )+1 ]
	input_rv_file = arguments[ arguments.index( '--rvin' )+1 ]

	output_fw_file = arguments[ arguments.index( '--fwout' )+1 ]
	output_rv_file = arguments[ arguments.index( '--rvout' )+1 ]
	
	
	#example: @NS500530:22:HVLCMBGX5:1:11101:10266:1041 1:N:0:CTTGTA
	#x: 1-99999
	#y: 1-99999
	#tile numer:1000-9999
	#lane: 1-8
	#flowcell ID: diverse string => get list of some
	#run number
	
	x = 1
	y = 1
	tile = 1000
	lane = 1
	flowcell = "HVLCMBGX5"
	run = 1
	
	with gzip.open( output_fw_file, "wb" ) as out_fw:
		with gzip.open( output_rv_file, "wb" ) as out_rv:
			with gzip.open( input_fw_file, "rb" ) as f1:
				with gzip.open( input_rv_file, "rb" ) as f2:
					line = f1.readline()
					while line:
						# --- generate read ID and update all variables --- #
						ID = ":".join( map( str, [ "@NS500530", run, flowcell, lane, tile, x, y ] ) )
						x += 1
						if x == 99999:
							x = 1
							y += 1
							if y == 99999:
								y = 1
								tile += 1
								print tile
								if tile == 9999:
									tile = 1000
									lane += 1
									print lane
									if lane == 8:
										lane = 1
										run += 1
										print run
										
						
						# --- write block in forward read file --- #
						out_fw.write( ID + ' 1:N:0:ACACAG\n' )
						out_fw.write( f1.readline() )
						f1.readline()
						out_fw.write( "+\n" )
						out_fw.write( f1.readline() )
						
						# --- write block in rv read file --- #
						f2.readline()
						out_rv.write( ID + ' 2:N:0:ACACAG\n' )
						out_rv.write( f2.readline() )
						f2.readline()
						out_rv.write( "+\n" )
						out_rv.write( f2.readline() )
	

if '--fwin' in sys.argv and '--rvin' in sys.argv and '--fwout' in sys.argv and '--rvout' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

