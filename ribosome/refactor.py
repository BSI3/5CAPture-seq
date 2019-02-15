import json
import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt


GENOME_SIZE = 23332839.0
MAX_CHROM_LENS = {"Pf3D7_01_v3": 640851,
	"Pf3D7_02_v3": 947102,
	"Pf3D7_03_v3": 1067971,
	"Pf3D7_04_v3": 1200490,
	"Pf3D7_05_v3": 1343557,
	"Pf3D7_06_v3": 1418242,
	"Pf3D7_07_v3": 1445207,
	"Pf3D7_08_v3": 1472805,
	"Pf3D7_09_v3": 1541735,
	"Pf3D7_10_v3": 1687656,
	"Pf3D7_11_v3": 2038340,
	"Pf3D7_12_v3": 2271494,
	"Pf3D7_13_v3": 2925236,
	"Pf3D7_14_v3": 3291936,
	"Pf3D7_API_v3": 34250,
	"Pf_M76611": 5967}
SELECTED_IDX_FILENAME = 'overlapped_extra_exonic'


def pre_process_depth(depth):
	for filename in depth:
		for chrom in depth[filename]:
			depth[filename][chrom]["genome_size"] = GENOME_SIZE

			for idx in range(1, MAX_CHROM_LENS[chrom] + 1):
				if idx not in depth[filename][chrom]["raw_depth"]:
					depth[filename][chrom]["raw_depth"][idx] = 0

	total_depth = 0
	for filename in depth:
		for chrom in depth[filename]:
			total_depth = total_depth + sum(depth[filename][chrom]["raw_depth"].values())

	for filename in depth:
		for chrom in depth[filename]:
			depth[filename][chrom]["total_depth"] = total_depth

	return depth

def read_depth(filename):
	depth_file = '/colossus/home/chadapohn/project_toner_reclu_compare/ribos_occ/raw_depth/' + filename + '.txt'

	depth = {}
	depth[filename] = {}

	for line in open(depth_file):
		chrom = line.split('\t')[0]
		idx = int(line.split('\t')[1])
		val = int(line.split('\t')[2].split('\n')[0])

		if chrom not in depth[filename]:
			depth[filename][chrom] = {}	
			depth[filename][chrom]["raw_depth"] = {}
			depth[filename][chrom]["raw_depth"][idx] = val

		if idx not in depth[filename][chrom]["raw_depth"]:
			depth[filename][chrom]["raw_depth"][idx] = val

	depth = pre_process_depth(depth)
	return depth

def normalize_depth(depth):
	for filename in depth:
		for chrom in depth[filename]:
			depth[filename][chrom]["normalize_depth"] = {}

			for idx in depth[filename][chrom]["raw_depth"]:
				depth[filename][chrom]["normalize_depth"][idx] = (depth[filename][chrom]["raw_depth"][idx] * depth[filename][chrom]["genome_size"]) / depth[filename][chrom]["total_depth"]

	return depth

# + "_" + SELECTED_IDX_FILENAME.split("_")[2]
def mark_occupancy(depth, filename):
	selected_idx_file = '/colossus/home/chadapohn/project_toner_reclu_compare/dinucleotide/work/compare_reclu_toner/position_weight_matrix/separate_out_in_gene/out_sep_4_groups/' + SELECTED_IDX_FILENAME + ".txt"
	log_file = "/colossus/home/chadapohn/project_toner_reclu_compare/ribos_occ/log_4_groups/" + SELECTED_IDX_FILENAME.split("_")[1] + "_" + SELECTED_IDX_FILENAME.split("_")[2] + "_" + filename + ".txt"
	log = open(log_file, "w")

	for line in open(selected_idx_file):
		if line[0] != '>':
				continue
		else:
			chrom = line.split(',')[0].split('>')[1]
			idx = int(line.split(',')[1])

			for filename in depth:	
				if "occupancy" not in depth[filename][chrom]:
					depth[filename][chrom]["occupancy"] = {}

				depth[filename][chrom]["occupancy"][idx] = {}
							
				lower_bound = idx - 1000
				upper_bound = idx + 1000

				for pos in range(lower_bound, upper_bound + 1): #iterate 2001 times
					if pos in depth[filename][chrom]["normalize_depth"]: 
						depth[filename][chrom]["occupancy"][idx][pos] = depth[filename][chrom]["normalize_depth"][pos]

					else:
						depth[filename][chrom]["occupancy"][idx][pos] = "x"

				
				log.write(chrom + "\t" + str(idx) + "\t" + str(depth[filename][chrom]["occupancy"][idx][idx]) + "\n") 
				
	return depth

def count_occupancy(depth):
	occupancy = {}
	for i in range(0, 2000 + 1): #iterate 2001 times
		occupancy[i] = []

	for filename in depth:
		for chrom in depth[filename]:
			if "occupancy" in depth[filename][chrom]:
				for idx in depth[filename][chrom]["occupancy"]:
					lower_bound = idx - 1000
					upper_bound = idx + 1000
					for i, pos in zip(range(0, 2000 + 1), range(lower_bound, upper_bound + 1)):
						occupancy[i].append(depth[filename][chrom]["occupancy"][idx][pos])

	y_oc = []
	y_sd = []

	for pos in occupancy:
		# print(pos)
		# print(occupancy[pos])
		total = []
		for idx in occupancy[pos]:
			if idx != "x":
				total.append(idx)
		
		avg = sum(total)/len(total)
		sd = np.std(total)

		# print(pos, sum(total), len(total), avg, sd)

		y_oc.append(avg)
		y_sd.append(sd)

	y_oc = np.array(y_oc)
	y_sd = np.array(y_sd)

	return(y_oc, y_sd)

# + ' ' + SELECTED_IDX_FILENAME.split("_")[2]
def plot_individual(depth, filename):
	y_oc, y_sd = count_occupancy(depth)

	x_oc = np.arange(-1000, 1001)
	y_oc = np.array(y_oc)
	y_sd = np.array(y_sd)

	plt.figure()
	plt.title('Average Ribosome Occupancy of ' + SELECTED_IDX_FILENAME.split("_")[1].upper() + ' ' + filename)
	plt.plot(x_oc, y_oc, 'blue', label = 'Average ' + filename)
	plt.plot(x_oc, y_oc + y_sd, 'gray', label = '+SD')
	plt.plot(x_oc, y_oc - y_sd, 'brown', label = '-SD')
	plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=3)
	plt.grid(True)
	plt.savefig('/colossus/home/chadapohn/project_toner_reclu_compare/ribos_occ/result/' + \
		SELECTED_IDX_FILENAME.split("_")[1] + '_' + filename + '.png', format ='png')
	plt.savefig('/colossus/home/chadapohn/project_toner_reclu_compare/ribos_occ/result/' + \
		SELECTED_IDX_FILENAME.split("_")[1] + '_' + filename + '.svg', format ='svg')
	plt.close()

#+ ' ' + SELECTED_IDX_FILENAME.split("_")[2].title()
# '_' + SELECTED_IDX_FILENAME.split("_")[2] +
def plot_overall(depths, filenames, colors):
	plt.figure()
	plt.title('Average Ribosome Occupancy of ' + SELECTED_IDX_FILENAME.split("_")[1].title() + ' ' + SELECTED_IDX_FILENAME.split("_")[2].title())

	for depth, filename, color in zip(depths, filenames, colors):
		y_oc, y_sd = count_occupancy(depth)

		x_oc = np.arange(-1000, 1001)
		y_oc = np.array(y_oc)
		y_sd = np.array(y_sd)

		plt.plot(x_oc, y_oc, color, label='Average ' + filename)

	plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True)
	ax = plt.gca()
	ax.set_ylim([0, 30])
	plt.grid(True)
	plt.savefig('/colossus/home/chadapohn/project_toner_reclu_compare/ribos_occ/result/' + \
		SELECTED_IDX_FILENAME.split("_")[1] + '_' + SELECTED_IDX_FILENAME.split("_")[2] + '.png', format ='png', bbox_inches="tight")
	plt.savefig('/colossus/home/chadapohn/project_toner_reclu_compare/ribos_occ/result/' + \
		SELECTED_IDX_FILENAME.split("_")[1] + '_' + SELECTED_IDX_FILENAME.split("_")[2] + '.svg', format ='svg', bbox_inches="tight")
	plt.close()


#'SRR2014726'
#'SRR2014727' 
#'SRR2014728'
#'SRR2014729'
#'SRR1378560'
#'SRR1378561'
#'SRR1378562'
#'SRR1378563'

filename_26 = 'SRR2014726'
depth_26 = read_depth(filename_26)
depth_26 = normalize_depth(depth_26)
depth_26 = mark_occupancy(depth_26, filename_26)
# for filename in depth_26:
# 	for chrom in depth_26[filename]:
# 		if "occupancy" in depth_26[filename][chrom]:
# 			for idx in depth_26[filename][chrom]["occupancy"]:
# 				print(chrom, idx)
			#print(depth_26[filename][chrom]["occupancy"][idx])

# for filename in depth_26:
# 	for chrom in depth_26[filename]:
# 		print(chrom, depth_26[filename][chrom]["total_depth"])
# plot_individual(depth_26, filename_26)

filename_27 = 'SRR2014727'
depth_27 = read_depth(filename_27)
depth_27 = normalize_depth(depth_27)
depth_27 = mark_occupancy(depth_27, filename_27)
#plot_individual(depth_27, filename_27)

filename_28 = 'SRR2014728'
depth_28 = read_depth(filename_28)
depth_28 = normalize_depth(depth_28)
depth_28 = mark_occupancy(depth_28, filename_28)
#plot_individual(depth_28, filename_28)

filename_29 = 'SRR2014729'
depth_29 = read_depth(filename_29)
depth_29 = normalize_depth(depth_29)
depth_29 = mark_occupancy(depth_29, filename_29)
#plot_individual(depth_29, filename_29)

filename_60 = 'SRR1378560'
depth_60 = read_depth(filename_60)
depth_60 = normalize_depth(depth_60)
depth_60 = mark_occupancy(depth_60, filename_60)
#plot_individual(depth_60, filename_60)

filename_61 = 'SRR1378561'
depth_61 = read_depth(filename_61)
depth_61 = normalize_depth(depth_61)
depth_61 = mark_occupancy(depth_61, filename_61)
#plot_individual(depth_61, filename_61)

filename_62 = 'SRR1378562'
depth_62 = read_depth(filename_62)
depth_62 = normalize_depth(depth_62)
depth_62 = mark_occupancy(depth_62, filename_62)
#plot_individual(depth_62, filename_62)

filename_63 = 'SRR1378563'
depth_63 = read_depth(filename_63)
depth_63 = normalize_depth(depth_63)
depth_63 = mark_occupancy(depth_63, filename_63)
#plot_individual(depth_63, filename_63)

depths = [depth_26, depth_27, depth_28, depth_29, depth_60, depth_61, depth_62, depth_63]
filenames = [filename_26, filename_27, filename_28, filename_29, filename_60, filename_61, filename_62, filename_63]
colors = ['#1B5E20', '#B71C1C', '#E65100', '#01579B', '#4CAF50', '#F44336', '#FF9800', '#03A9F4']
plot_overall(depths, filenames, colors)