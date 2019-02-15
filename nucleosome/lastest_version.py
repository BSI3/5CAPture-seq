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
SELECTED_IDX_FILENAME = 'overlapped_cec2'


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
	depth_file = '/colossus/home/chadapohn/project_toner_reclu_compare/nucl_occ/work/bwa/mapped_read/raw_depth/' + filename + '.txt'

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

def normalize_control(depth, control, filename):
	log_file = "/colossus/home/chadapohn/project_toner_reclu_compare/nucl_occ/work/bwa/mapped_read/normalized_depth/" + filename + ".txt"
	log = open(log_file, "w")

	for filename in depth:
		for chrom in depth[filename]:
			for idx in depth[filename][chrom]["normalize_depth"]:
				norm = (depth[filename][chrom]["normalize_depth"][idx] + 0.001) / (control["control"][chrom]["normalize_depth"][idx] + 0.001)

				log.write(chrom + "\t" + str(idx) + "\t" + str(norm) + "\n") 


def mark_occupancy(depth, control, filename):
	selected_idx_file = '/colossus/home/chadapohn/project_toner_reclu_compare/dinucleotide/work/compare_reclu_toner/position_weight_matrix/separate_out_in_gene/CEC_nuc_ribo_ChIP_1000starts/' + SELECTED_IDX_FILENAME + '.txt'

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
						if (depth[filename][chrom]["normalize_depth"][pos] == 0) & (control["control"][chrom]["normalize_depth"][pos] == 0):
							depth[filename][chrom]["occupancy"][idx][pos] = "x"
							# print(chrom, idx, pos, depth[filename][chrom]["normalize_depth"][pos], control["control"][chrom]["normalize_depth"][pos], depth[filename][chrom]["occupancy"][idx][pos])

						else:
							depth[filename][chrom]["occupancy"][idx][pos] = (depth[filename][chrom]["normalize_depth"][pos] + 0.001) / (control["control"][chrom]["normalize_depth"][pos] + 0.001)
							# print(chrom, idx, pos, depth[filename][chrom]["normalize_depth"][pos] + 0.001, control["control"][chrom]["normalize_depth"][pos] + 0.001, depth[filename][chrom]["occupancy"][idx][pos])

					else:
						depth[filename][chrom]["occupancy"][idx][pos] = "x"
						# print(chrom, idx, pos, depth[filename][chrom]["occupancy"][idx][pos])
		
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

def plot_individual(depth, filename):
	y_oc, y_sd = count_occupancy(depth)

	x_oc = np.arange(-1000, 1001)
	y_oc = np.array(y_oc)
	y_sd = np.array(y_sd)

	plt.figure()
	plt.title('Average Ribosome Occupancy of ' + SELECTED_IDX_FILENAME.split("_")[2].title() + ' ' + SELECTED_IDX_FILENAME.split("_")[1] + ' ' + filename)

	plt.plot(x_oc, y_oc, 'blue', label = 'Average ' + filename)
	plt.plot(x_oc, y_oc + y_sd, 'gray', label = '+SD')
	plt.plot(x_oc, y_oc - y_sd, 'brown', label = '-SD')
	plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=3)
	plt.grid(True)
	plt.savefig('/colossus/home/chadapohn/project_toner_reclu_compare/nucl_occ/work/bwa/mapped_read/result/' + \
		SELECTED_IDX_FILENAME.split("_")[2] + '_' + SELECTED_IDX_FILENAME.split("_")[1] + '_' + filename + '.png', format ='png')
	plt.savefig('/colossus/home/chadapohn/project_toner_reclu_compare/nucl_occ/work/bwa/mapped_read/result/' + \
		SELECTED_IDX_FILENAME.split("_")[2] + '_' + SELECTED_IDX_FILENAME.split("_")[1] + '_' + filename + '.svg', format ='svg')
	plt.close()


def plot_overall(depths, filenames, colors):
	plt.figure()
	plt.title('Average Nucleosome Occupancy of ' + SELECTED_IDX_FILENAME.split("_")[1].upper())

	for depth, filename, color in zip(depths, filenames, colors):
		y_oc, y_sd = count_occupancy(depth)

		x_oc = np.arange(-1000, 1001)
		y_oc = np.array(y_oc)
		y_sd = np.array(y_sd)

		plt.plot(x_oc, y_oc, color, label='Average ' + filename.title())

	plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True)
	ax = plt.gca()
	ax.set_ylim([0.8, 1.2])
	plt.grid(True)
	plt.savefig('/colossus/home/chadapohn/project_toner_reclu_compare/nucl_occ/work/bwa/mapped_read/result/' + \
		SELECTED_IDX_FILENAME.split("_")[1] + '.png', format ='png', bbox_inches="tight")
	plt.savefig('/colossus/home/chadapohn/project_toner_reclu_compare/nucl_occ/work/bwa/mapped_read/result/' + \
		SELECTED_IDX_FILENAME.split("_")[1] + '.svg', format ='svg', bbox_inches="tight")
	plt.close()


filename_control = 'control'
depth_control = read_depth(filename_control)
depth_control = normalize_depth(depth_control)

filename_t5 = 't5'
depth_t5 = read_depth(filename_t5)
depth_t5 = normalize_depth(depth_t5)
depth_t5 = normalize_control(depth_t5, depth_control, filename_t5)
# depth_t5 = mark_occupancy(depth_t5, depth_control, filename_t5)

filename_t10 = 't10'
depth_t10 = read_depth(filename_t10)
depth_t10 = normalize_depth(depth_t10)
depth_t10 = normalize_control(depth_t10, depth_control, filename_t10)
# depth_t10 = mark_occupancy(depth_t10, depth_control, filename_t10)

filename_t15 = 't15'
depth_t15 = read_depth(filename_t15)
depth_t15 = normalize_depth(depth_t15)
depth_t15 = normalize_control(depth_t15, depth_control, filename_t15)
# depth_t15 = mark_occupancy(depth_t15, depth_control, filename_t15)

filename_t20 = 't20'
depth_t20 = read_depth(filename_t20)
depth_t20 = normalize_depth(depth_t20)
depth_t20 = normalize_control(depth_t20, depth_control, filename_t20)
# depth_t20 = mark_occupancy(depth_t20, depth_control, filename_t20)

filename_t25 = 't25'
depth_t25 = read_depth(filename_t25)
depth_t25 = normalize_depth(depth_t25)
depth_t25 = normalize_control(depth_t25, depth_control, filename_t25)
# depth_t25 = mark_occupancy(depth_t25, depth_control, filename_t25)

filename_t30 = 't30'
depth_t30 = read_depth(filename_t30)
depth_t30 = normalize_depth(depth_t30)
depth_t30 = normalize_control(depth_t30, depth_control, filename_t30)
# depth_t30 = mark_occupancy(depth_t30, depth_control, filename_t30)

filename_t35 = 't35'
depth_t35 = read_depth(filename_t35)
depth_t35 = normalize_depth(depth_t35)
depth_t35 = normalize_control(depth_t35, depth_control, filename_t35)
# depth_t35 = mark_occupancy(depth_t35, depth_control, filename_t35)

filename_t40 = 't40'
depth_t40 = read_depth(filename_t40)
depth_t40 = normalize_depth(depth_t40)
depth_t40 = normalize_control(depth_t40, depth_control, filename_t40)
# depth_t40 = mark_occupancy(depth_t40, depth_control, filename_t40)

# depths = [depth_t5, depth_t10, depth_t15, depth_20, depth_25, depth_30, depth_35, depth_40]
# filenames = [filename_t5, filename_t10, filename_t15, filename_t20, filename_t25, filename_t30, filename_t35, filename_t40]
# colors = ['#C62828', '#2E7D32', '#304FFE', '#EF6C00', '#AA00FF', '#FFAB00', '#795548', '#9E9E9E']
# plot_overall(depths, filenames, colors)


