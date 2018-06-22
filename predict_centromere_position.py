import matplotlib.pyplot as plt
import heapq
import sys, os

__usage__ = """ python centromere_prediction.py
			--input_vcf <INPUT_VCF_FILE,DIRECTORY_FOR_OUTPUT>\n
			--input_fasta <INPUT_REFERENCE_FILE>\n
			--centromere_pos  <OPTIONAL,ALL_CENTROMERE_POSITIONS_FROM_CHR_1_TO_N_SEPERATED_BY_COMMA,SEVERAL_POSSIBLE_CENTROMERES_FOR_ONE_CHROMOSOME_CAN_BE_SEPARATED_BY_PLUS>\n
			\n
			Example for --centromere_pos with four chromosomes: \n
			--centromere_pos 13400000,12000000+11500000,14000000,12500000 \n
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

def load_variants_from_vcf(vcf_file, number_of_chr, chr_lengths, resolution_factor):
	"""! @brief loads the variant information from a VCF file """
	snps_per_chr = []
	indels_per_chr = []
	qualities = []
	for i in range(number_of_chr):
		snps_per_chr.append([])
	for i in range(number_of_chr):
		indels_per_chr.append([])
	for i in range(number_of_chr):
		qualities.append([])
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
					#append positions and quality
					quality = (int(parts[1]), float(parts[5]))
					qualities[chr_nr - 1].append(quality)
				except:
					line = f.readline()
					pass	
			line = f.readline()
	snp_data = []
	for i in range(number_of_chr):
		snp_data.append([])
	indel_data = []
	for i in range(number_of_chr):
		indel_data.append([])
	quality_data = []
	for i in range(number_of_chr):
		quality_data.append([])
	for idx, chr_length in enumerate(chr_lengths):
		#identify positions of variants
		resolution = (chr_length / float(100)) * resolution_factor
		upper_lim = resolution
		lower_lim = 0
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

				quality_sum = 0
				quality_number = 0
				for (pos, qual) in qualities[ idx ]:
					if pos <= upper_lim and pos > lower_lim:
						quality_sum = quality_sum + qual
						quality_number = quality_number + 1
				try:
					quality_data[idx].append(quality_sum / quality_number)
				except:
					quality_data[idx].append(0)
					pass

				if len(snp_tmp) == 0:
					snp_data[idx].append(1)		#to avoid dividing by 0
				else:
					snp_data[idx].append(len(snp_tmp))
				if len(indel_tmp) == 0:
					indel_data[idx].append(1)	#to avoid dividing by 0
				else:
					indel_data[idx].append(len(indel_tmp))
			upper_lim += resolution
			lower_lim += resolution
	return snp_data, indel_data, quality_data

def calculate_relative_SNP_to_inDel_ratio(number_of_chr, chr_lengths, snp_data, indel_data, quality_data, resolution_factor):
	"""! @brief calculate the relative SNP to InDel ratio and prepare plotting """
	#fill perc_snps and perc_indels
	perc_snps = []	#percent values from snp_data
	for i in range(number_of_chr):
		perc_snps.append([])
	perc_indels = [] #percent values from indel_data
	for i in range(number_of_chr):
		perc_indels.append([])	
	for i, chrom in enumerate(snp_data):
		indiv_max = max(chrom)
		temp = []
		for j, pos in enumerate(chrom):
			temp.append(float(pos) / float(indiv_max))
		perc_snps[i] = temp
	for i, chrom in enumerate(indel_data):
		indiv_max = max(chrom)
		temp = []
		for j, pos in enumerate(chrom):
			 temp.append(float(pos) / float(indiv_max))
		perc_indels[i] = temp		
	
	ratios = []
	for i in range(number_of_chr):
		ratios.append([])
	xdots = []
	for i in range(number_of_chr):
		xdots.append([])
	
	for c in range(0,number_of_chr):
		ratio = []
		for i in range(0,len(perc_snps[c])):
			r = perc_snps[c][i] / perc_indels[c][i]
			r = round(r,2)
			ratio.append(r)
		ratios[c] = ratio
		
		#define x intervals for plotting
		x_axis = []
		resolution = (chr_lengths[c] / float(100)) * resolution_factor
		steps = float(resolution) / 1000000
		q = steps
		for x in range(0,len(perc_snps[c])):
			x_axis.append(q)
			q = q + steps	
		xdots[c] = x_axis
	plots = zip(xdots, ratios, quality_data) 	#zip the chromosomes together
	return plots

def loop_plot(plots, fig_prefix, chr_lengths, centromere_pos):
	"""! @brief for each chromosome plot relative SNP/InDel ratio and find centromere position"""
	figs = {}
	axs = {}
	qual = {}
	all_fig_names = ""
	for idx, plot in enumerate(plots):
		figs[idx] = plt.figure(figsize = (10,5))
		axs[idx] = plt.subplot()
		qual[idx] = axs[idx].twinx()
		
		axs[idx].set_xlim(right = chr_lengths[idx] / float(1000000), left=0, auto = False)
		a, = axs[idx].plot(plot[0], plot[1], label = "Relative SNP/InDel\nratio", color = "blue", linewidth = 2)
		q,= qual[idx].plot(plot[0],plot[2], color = "red", label = "Quality")

		#get highest peak
		high = heapq.nlargest(2, plot[1])
		index_one = plot[1].index(high[0])
		index_two = plot[1].index(high[1])
		cen_one = plot[0][index_one]
		cen_two = plot[0][index_two]
		chr_length = chr_lengths[idx] / float(1000000)
		if chr_length - cen_one < 1 or chr_length - cen_one > chr_length - 1: #avoid peaks 1 Mbp near both ends of the chromosome
			i = cen_two
		else:
			i = cen_one
		
		ymin, ymax = qual[idx].get_ylim()
		xmin, xmax = qual[idx].get_xlim()
		c, = qual[idx].plot((i, i), (ymin, ymax), linestyle = "--", color = "xkcd:kelly green", zorder = 3, label = "Predicted\ncentromere pos.")
		qual[idx].text(i, ymax, "{0:.1f}".format(i), ha = "center", va = "bottom")

		#append ticks and labels every 5 Mbp
		ticks = []
		labels = []
		curr_mb = 0
		while True:
			if curr_mb > xmax:
				break
			else:
				ticks.append(curr_mb)
				labels.append(curr_mb)
				curr_mb = curr_mb + 5
		axs[idx].set_xticks(ticks)
		axs[idx].set_xticklabels(labels)
		
		#if given, append reference centromere positions
		if centromere_pos is not None:
			for p in centromere_pos[idx]:					
					d, = qual[idx].plot((p, p),(ymin, ymax), "-", color = "black", label = "Reference\ncentromere pos.")

		axs[idx].plot((xmin, xmax), (1, 1), color = "orange") 
		axs[idx].set_title("Relative SNP/InDel ratio over chr" + str(idx + 1))
		axs[idx].set_xlabel("[ Mbp ]")
		axs[idx].set_ylabel("Relative SNP/InDel ratio")
		qual[idx].set_ylabel("Quality")
		qual[idx].yaxis.label.set_color(q.get_color())
		qual[idx].tick_params(axis = "y", colors = q.get_color())
		axs[idx].yaxis.label.set_color(a.get_color())
		axs[idx].tick_params(axis = "y", colors = a.get_color())

		# Shrink current axis by 20%
		box = axs[idx].get_position()
		axs[idx].set_position([box.x0, box.y0, box.width * 0.8, box.height])
		qual[idx].set_position([box.x0, box.y0, box.width * 0.8, box.height])

		if centromere_pos is not None:
			axs[idx].legend(handles = [a,c,d,q],loc='center left', bbox_to_anchor=(1.1, 0.5))
		else:
			axs[idx].legend(handles = [a,c,q],loc='center left', bbox_to_anchor=(1.1, 0.5))

		figs[idx].savefig(fig_prefix + "_chr" + str(idx + 1) + ".png", dpi=300)
		all_fig_names = all_fig_names + fig_prefix + "_chr" + str(idx + 1) + ".png" + " "

	cmd = "montage " + all_fig_names + " -geometry +2+2 -tile 2x " + fig_prefix + "_all_chr.png"
	os.popen(cmd)
	return figs, axs

def main(arguments):
	"""! @brief gets all arguments from command line and calls functions for computing and plotting SNP/InDel ratio """
	resolution_factor = 2
	vcf_file = arguments[arguments.index("--input_vcf" ) + 1]
	fasta_file = arguments[arguments.index("--input_fasta" ) + 1]
	directory = vcf_file.strip().split("/")
	del directory[-1]
	fig_prefix = "/".join(directory) + "/centromere_prediction"
	number_of_chr, chr_lengths = load_chr_number_and_lengths_from_fasta(fasta_file)
	print "done reading fasta"
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
	
	snp_data, indel_data, quality_data = load_variants_from_vcf(vcf_file, number_of_chr, chr_lengths, resolution_factor)
	print "done reading vcf file"
	plots = calculate_relative_SNP_to_inDel_ratio(number_of_chr, chr_lengths, snp_data, indel_data, quality_data, resolution_factor)
	figs, axs = loop_plot(plots, fig_prefix, chr_lengths, centromere_pos)
	plt.show()
	
if __name__ == "__main__":

	if "--input_vcf" in sys.argv and "--input_fasta" in sys.argv:
		main(sys.argv)
	else:
		sys.exit(__usage__)
