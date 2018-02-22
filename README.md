# script_collection
Collection of scripts to solve small bioinformatic challenges.


identify_RBHs.py
Identification of Reciprocal Best BLAST Hits (RBHs) between to sets of sequences (protein/DNA). The scripts constructs BLAST databases and runs blastp/blastn in both directions. RBHs are identified and writen to a text file ('RBH_file.txt') in the specified output directory.
Requirements:
1) Python 2.7.x (other Python 2 versions should work as well)
2) BLAST (makeblastdb, blastn, and blastp should be in PATH)
Usage:
python identify_RBHs.py \
--input1 <FASTA_FILE_1> \
--input2 <FASTA_FILE2> \
--prefix <OUTPUT_DIRECTORY_NAME> \
--seq_type <prot|nucl> \
--cpu <INT>
