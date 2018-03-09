# script_collection
Collection of scripts to solve small bioinformatic challenges.


#identify_RBHs.py

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
--cpu <NUMBER_OF_CPUs_TO_USE>


Suggested citation:

Pucker et al., 2016: 'A De Novo Genome Sequence Assembly of the Arabidopsis thaliana Accession Niederzenz-1 Displays Presence/Absence Variation and Strong Synteny'
http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0164321





#sort_contigs_on_ref.py

Whole Genome Shotgun (WGS) assembly contigs can be ordered and oriented based on an available reference sequence. This script does a placement of all given sequences based on the central position of their best BLASTn hit against the reference sequence. A new FASTA file is constructed, in which all seqeuences are saved under new systematic names (scaffold<running_number>). Association between old and new names is printed during this process and can easily be written into a documentation file.

Requirements:

1) Python 2.7.x (other Python 2 versions should work as well)
2) BLAST (makeblastdb and blastn should be in PATH)

Usage:

python sort_contigs_on_ref.py \
--contig_file <FULL_PATH_TO_FILE> \
--ref_file <FULL_PATH_TO_FILE> \
--output_dir <FULL_PATH_TO_DIR> > <DOCUMENTATION_FILE>

Suggested citation:

Pucker et al., 2016: 'A De Novo Genome Sequence Assembly of the Arabidopsis thaliana Accession Niederzenz-1 Displays Presence/Absence Variation and Strong Synteny'
http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0164321





#split_FASTQ.py

Splits FASTQ file with alternating mate1 and mate2 reads of paired-end sequencing into two separate files with mate1 and mate2, respectively. Can be applied after downloading FASTQ files from the SRA via webbrowser. New files will be placed next to the original file with '_1' and '_2' added to their file base name. This script can handle raw FASTQ files (.fastq) as well as gzip compressed files (.fastq.gz).

Requirements:

1) Python 2.7.x (other Python 2 versions should work as well)

Usage:

python split_FASTQ.py \
--in_file <FULL_PATH_TO_FILE>


Suggested citation:

this repository





#sort_vcf_by_fasta.py

A given VCF file is sorted based on the provided FASTA file. The chromosome order and numeric positions within the chromosome sequences are taken into account to adjust the VCF file. This can be helpful during variant calling with GATK.

Requirements:

1) Python 2.7.x (other Python 2 versions should work as well)

Usage:

python sort_vcf_by_fasta.py \
--vcf <FULL_PATH_TO_INPUT_VCF> \
--fasta <FULL_PATH_TO_INPUT_FASTA_FILE> \
--output <FULL_PATH_TO_OUTPUT_VCF_FILE>


Suggested citation:

this repository








