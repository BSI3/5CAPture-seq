import json
import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt


MAX_CHROM_LENS = {"M76611": 5967,
"PFC10_API_IRAB": 34242,
"Pf3D7_01_v3": 640851,
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
"Pf3D7_14_v3": 3291936}
SELECTED_IDX_FILENAME = 'overlapped_cec2'


def pre_process_depth(depth):
	for filename in depth:
		for chrom in depth[filename]:

			for idx in range(1, MAX_CHROM_LENS[chrom] + 1):
				if idx not in depth[filename][chrom]["raw_depth"]:
					depth[filename][chrom]["raw_depth"][idx] = 0

	return depth

def read_depth(filename):
	depth_file = '/colossus/home/chadapohn/project_toner_reclu_compare/rna_polymerase/ChIPSeq/' + filename + '.txt'

	depth = {}
	depth[filename] = {}

	for line in open(depth_file):
		chrom = line.split('\t')[0]
		idx = int(line.split('\t')[1])
		val = float(line.split('\t')[2].split('\n')[0])

		if chrom not in depth[filename]:
			depth[filename][chrom] = {}	
			depth[filename][chrom]["raw_depth"] = {}
			depth[filename][chrom]["raw_depth"][idx] = val

		if idx not in depth[filename][chrom]["raw_depth"]:
			depth[filename][chrom]["raw_depth"][idx] = val

	depth = pre_process_depth(depth)
	return depth

def mark_occupancy(depth, filename):
	selected_idx_file = '/colossus/home/chadapohn/project_toner_reclu_compare/dinucleotide/work/compare_reclu_toner/position_weight_matrix/separate_out_in_gene/CEC_nuc_ribo_ChIP_1000starts/' + SELECTED_IDX_FILENAME + ".txt"
	#log_file = "/colossus/home/chadapohn/project_toner_reclu_compare/rna_polymerase/log_CEC_nuc_ribo_ChIP_1000starts/" + SELECTED_IDX_FILENAME.split("_")[1] + "_" + filename + ".txt"
	#log = open(log_file, "w")

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
					if pos in depth[filename][chrom]["raw_depth"]: 
						depth[filename][chrom]["occupancy"][idx][pos] = depth[filename][chrom]["raw_depth"][pos]

					else:
						depth[filename][chrom]["occupancy"][idx][pos] = "x"

				
				#log.write(chrom + "\t" + str(idx) + "\t" + str(depth[filename][chrom]["raw_depth"][idx][idx]) + "\n") 
				
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
		total = []
		for idx in occupancy[pos]:
			if idx != "x":
				total.append(idx)
		
		avg = sum(total)/len(total)
		sd = np.std(total)

		y_oc.append(avg)
		y_sd.append(sd)

	y_oc = np.array(y_oc)
	y_sd = np.array(y_sd)

	return(y_oc, y_sd)

def plot_overall(depths, filenames, colors):
	plt.figure()
	plt.title('Average RNApolII Occupancy of ' + SELECTED_IDX_FILENAME.split("_")[1].upper())

	for depth, filename, color in zip(depths, filenames, colors):
		y_oc, y_sd = count_occupancy(depth)

		x_oc = np.arange(-1000, 1001)
		y_oc = np.array(y_oc)
		y_sd = np.array(y_sd)

		plt.plot(x_oc, y_oc, color, label='Average ' + filename)

	plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True)
	ax = plt.gca()
	ax.set_ylim([0, 14])
	plt.grid(True)
	plt.savefig('/colossus/home/chadapohn/project_toner_reclu_compare/rna_polymerase/result/' + \
		SELECTED_IDX_FILENAME.split("_")[1] + '.png', format ='png', bbox_inches="tight")
	plt.savefig('/colossus/home/chadapohn/project_toner_reclu_compare/rna_polymerase/result/' + \
		SELECTED_IDX_FILENAME.split("_")[1] + '.svg', format ='svg', bbox_inches="tight")
	plt.close()


filename_er = 'GSE85478_ChIPSeq_ER_normalized'
depth_er = read_depth(filename_er)
depth_er = mark_occupancy(depth_er, filename_er)

filename_et = 'GSE85478_ChIPSeq_ET_normalized'
depth_et = read_depth(filename_et)
depth_et = mark_occupancy(depth_et, filename_et)

filename_ls = 'GSE85478_ChIPSeq_LS_normalized'
depth_ls = read_depth(filename_ls)
depth_ls = mark_occupancy(depth_ls, filename_ls)


depths = [depth_er, depth_et, depth_ls]
filenames = [filename_er, filename_et, filename_ls]
colors = ['#2E7D32', '#C62828', '#304FFE']
plot_overall(depths, filenames, colors)