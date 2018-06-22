import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys, math

__usage__ = """ python genome_wide_distribution_of_variants.py
			--input_vcf <INPUT_VCF_FILE,DIRECTORY_FOR_OUTPUT>\n
			--input_fasta <INPUT_REFERENCE_FILE>\n
			--resolution <INT_FOR_WINDOW_SIZE>\n
			--centromere_pos  <OPTIONAL,ALL_CENTROMERE_POSITIONS_FROM_CHR_1_TO_N_SEPERATED_BY_COMMA,SEVERAL_POSSIBLE_CENTROMERES_FOR_ONE_CHROMOSOME_CAN_BE_SEPARATED_BY_PLUS>\n
			\n
			Example for --centromere_pos with four chromosomes: \n
			--centromere_pos 13400000,12000000+11500000,14000000,12500000 \n
			Resolution 100000 is recommended.\n
			"""

def load_chr_number_and_lengths_from_fasta(fasta_file):
	"""! @brief load the number of chromosomes and all chromosome lengths from FASTA file """
	chr_lengths = []
	with open(fasta_file, "r" ) as f:
		line = f.readline()
		length = 0
		chr_nr = ""
		while line:
			if line[0] == ">":
				parts = line.strip().split(" ")
				if "chr" in line:
					chr_nr = parts[0].replace(">chr", "")
				else:
					chr_nr = parts[0].replace(">Chr", "")
				try:
					chr_nr = int(chr_nr)
					length = 0
					line = f.readline()
					while line:
						if line[0] == ">":
							chr_lengths.append((chr_nr, length))
							break
						length = length + len(line)
						line = f.readline()
				except:
					line = f.readline()
					pass
			else:
				line = f.readline()
	chr_lengths.sort(key = lambda s: s[0])
	for i, elem in enumerate(chr_lengths):
		chr_lengths[i] = elem[1]
	number_of_chr = len(chr_lengths)
	return number_of_chr, chr_lengths

def load_variants_from_vcf(vcf_file, number_of_chr):
	"""! @brief loads the variant information from a VCF file """
	snps_per_chr = []
	indels_per_chr = []
	for i in range(number_of_chr):
		snps_per_chr.append([])
	for i in range(number_of_chr):
		indels_per_chr.append([])
	with open(vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != "#":
				try:
					parts = line.strip().split("\t")
					name = str(parts[0])
					if "chr" in name:
						name = name.replace("chr", "")
					elif "Chr" in name:
						name = name.replace("Chr", "")
					chr_nr = int(name)
					if len(parts[3]) == len(parts[4]):	#SNP 
						snps_per_chr[chr_nr - 1].append(int(parts[1]))
					else:					#InDel
						indels_per_chr[chr_nr - 1].append(int(parts[1]))
				except:
					line = f.readline()
					pass	
			line = f.readline()
	return snps_per_chr, indels_per_chr

def construct_plot(snps_per_chr, indels_per_chr, number_of_chr, chr_lengths, result_file, result_table, resolution, centromere_pos):
	"""! @brief construct variant genome distribution plot """
	
	f = plt.figure(figsize = (15,15))
	
	#if more than 10 chromosomes should be displayed, make two columns
	if number_of_chr >= 10:	
		gs = gridspec.GridSpec(int(math.ceil(float(number_of_chr) / 2)), 2)	#make two columns with each half of the chromosomes
	else:
		gs = gridspec.GridSpec(number_of_chr, 1)

	#length of the x-axis is the individual chromosome length in Mbp
	x_lim = []
	for index, chr_len in enumerate(chr_lengths):
		x_lim.append(int(str(chr_len)[:2]))

	with open(result_table, "w") as out:
		for idx, chr_length in enumerate(chr_lengths):
			chr_name = "Chr" + str(idx + 1)
			
			#identify positions of variants
			upper_lim = resolution
			lower_lim = 0
			
			snp_data = []
			indel_data = []
			while True:
				if upper_lim >= chr_length:
					break
				else:
					snp_tmp = []
					indel_tmp = []
					for SNP in snps_per_chr[idx]:
						if SNP <= upper_lim and SNP > lower_lim:
							snp_tmp.append("X")
					for indel in indels_per_chr[ idx ]:
						if indel <= upper_lim and indel > lower_lim:
							indel_tmp.append("X")
					snp_data.append(len(snp_tmp))
					indel_data.append(len(indel_tmp))
				upper_lim += resolution
				lower_lim += resolution
			
			#length of the left y-axis is the max. number of SNPs		
			y_max_lim = max(snp_data)
			if max(indel_data) > y_max_lim:
				y_max_lim = max(indel_data)

			#get max of x-axis
			display_x_max = max(x_lim)

			#placing the plots and creating two y-axis
			if number_of_chr >= 10:	#plot in two columns
				if idx >= int(math.ceil(float(number_of_chr) / 2)):	#plot in the second column
					#for all chromosomes the x-axis length is relative to length of longest chromosome
					current_grid = gridspec.GridSpecFromSubplotSpec(1, display_x_max, subplot_spec = gs[idx-int(math.ceil(float(number_of_chr) / 2)),1])
					ax_b = plt.subplot(current_grid[0: x_lim[idx]])
					ax_a = ax_b.twinx()
				else:	#plot in the first column
					current_grid = gridspec.GridSpecFromSubplotSpec(1, display_x_max, subplot_spec = gs[idx,0])
					ax_b = plt.subplot(current_grid[ 0: x_lim[idx]])
					ax_a = ax_b.twinx()
			else:	#plot in one column
				#for all chromosomes the x-axis length is relative to length of longest chromosome
				current_grid = gridspec.GridSpecFromSubplotSpec(1, display_x_max, subplot_spec = gs[idx,0])
				ax_b = plt.subplot(current_grid[ 0: x_lim[idx]])
				ax_a = ax_b.twinx()

			#plotting SNP and InDel distribution
			res = resolution / float(1000000)
			start = resolution / float(1000000)
			ratios = []
			x = []
			max_snp = max(snp_data)
			max_indel = max(indel_data)
			for i, snps in enumerate(snp_data):
				a, = ax_a.plot(((res, res)), (0, snps), "-", color = "black", label = "SNPs")
				b, = ax_b.plot(((res + 0.1, res + 0.1)), (0, indel_data[i]), "-", color = "xkcd:dull red", label = "InDels")
				res = res + start
			res = resolution / float(1000000)
			start = resolution / float(1000000)
			ax_indel = ax_b.twinx()
			for i, snps in enumerate(snp_data):
				if (float(snps)/max_snp) > (float(indel_data[i])/max_indel):
					h, = ax_indel.plot(((res + 0.1, res + 0.1)), (0, indel_data[i]), "-", color = "xkcd:dull red")
				res = res + start

			#plotting centromere positions if given
			if centromere_pos is not None:
				for p in centromere_pos[idx]:					
					c, = ax_indel.plot((p, p),(0, y_max_lim), "-", color = "xkcd:kelly green", linewidth = 3, zorder = 3, label = "Centromere")

			#writing data into output table
			out.write("Chr" + str(idx + 1) + "SNPs:\t" + "\t".join(map(str, snp_data)) + "\n")
			out.write("Chr" + str(idx + 1) + "InDels:\t" + "\t".join(map(str, indel_data)) + "\n")
			
			#adding and correcting labels, legend and size
			current_labels = ax_a.get_xticks()
			labels = []
			for each in current_labels:
				labels.append(each)
			ax_a.set_xticks([])
			ax_a.set_xticks(labels)
			ax_a.set_title(chr_name)
			ax_b.set_xlabel("[ Mbp ]")
			ax_a.set_ylabel("Number of SNPs")
			ax_a.axis([0, x_lim[idx], 0, y_max_lim])		
			ax_b.tick_params("y", colors = "xkcd:dull red")
			ax_b.set_ylabel("Number of InDels", color = "xkcd:dull red")
			ax_a.set_xlim(right = x_lim[idx], auto = False)
			ax_b.set_ylim(bottom = 0, auto = False)	
			ymin, ymax = ax_b.get_ylim()
			ax_indel.set_ylim(bottom = 0, top = ymax, auto = False)
			ax_indel.axes.get_yaxis().set_visible(False)

	plt.figlegend([a], ["SNPs"],(0.20, 0.0), prop = {"size": "x-large"}, frameon = False)
	plt.figlegend([b], ["Indels"], (0.30, 0.0), prop = {"size": "x-large"}, frameon = False)
	if centromere_pos is not None:
		plt.figlegend([c], ["Centromere"], (0.40, 0.0), prop = {"size": "x-large"}, frameon = False)

	gs.tight_layout(f)	
	gs.update(hspace = 0.8, wspace = 0.4)
	f.set_figheight(40)
	f.set_figwidth(40)
	plt.show()
	f.savefig(result_file, dpi = 300)

def main(arguments):
	"""! @brief gets all arguments from command line and calls functions for computing and plotting variant distribution """
	resolution = float(arguments[arguments.index("--resolution") + 1])
	vcf_file = arguments[arguments.index("--input_vcf" ) + 1]
	fasta_file = arguments[arguments.index("--input_fasta") + 1]
	directory = vcf_file.strip().split("/")
	del directory[-1]
	result_file = "/".join(directory) + "/genome_wide_small_variants_" + str(int(resolution)) + ".png"
	result_table = "/".join(directory) + "/genome_wide_small_variants_" + str(int(resolution)) + ".txt"
	number_of_chr, chr_lengths = load_chr_number_and_lengths_from_fasta(fasta_file)
	if "--centromere_pos" in arguments:
		centromere_pos = []
		centromere_pos_raw = arguments[arguments.index("--centromere_pos") + 1]
		for i in range(number_of_chr):
			centromere_pos.append([])
		chromosomes = centromere_pos_raw.split(",")
		for i, chromosome in enumerate(chromosomes):
			centromere = chromosome.split("+")
			for idx, each in enumerate(centromere):
				centromere[idx] = float(each) / float(1000000)
			centromere_pos[i] = centromere
	else:
		centromere_pos = None
	snps_per_chr, indels_per_chr = load_variants_from_vcf(vcf_file, number_of_chr)
	construct_plot(snps_per_chr, indels_per_chr, number_of_chr, chr_lengths, result_file, result_table, resolution, centromere_pos)

if __name__ == "__main__":

	if "--input_vcf" in sys.argv and "--input_fasta" in sys.argv and "--resolution" in sys.argv:
		main(sys.argv)
	else:
		sys.exit(__usage__)
