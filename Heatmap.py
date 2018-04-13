import sys
import seaborn as sns
import argparse
import numpy as np
import scipy.spatial
from matplotlib import pyplot as plt

def readData(expfile):
	IF=open(expfile,'r')
	line=IF.readline()
	lst_sample=line.strip().split('\t')[1:]
			
	lst_gene=[]
	lst_exps=[]
	for line in IF:
		s=line.strip().split('\t')
		gene, lst_exp = s[0], map(float,s[1:])
		lst_gene.append(gene)
		lst_exps.append(lst_exp)
	arr_exp = np.array(lst_exps)
	return lst_sample, lst_gene, arr_exp

def label2struct(lst_sample):
	dic_st2i_val={}
	for i, sample in enumerate(lst_sample):
		if ',' in sample:
			sample = sample.split(',')[1]
		s_val, t_val = sample.split('_')
		s_val, t_val = int(s_val.strip('s')), int(t_val.strip('t'))
		if s_val not in dic_st2i_val:
			dic_st2i_val[s_val]={}
		if t_val not in dic_st2i_val[s_val]:
			dic_st2i_val[s_val][t_val]=i

	lst_st2idx=[]
	for s, s_val in enumerate(sorted(dic_st2i_val.keys(),key=lambda x:float(x))):
		for t, t_val in enumerate(sorted(dic_st2i_val[s_val].keys(),key=lambda x:float(x))):
			if len(lst_st2idx) <= s:
				lst_st2idx.append([])
			if len(lst_st2idx[s]) <= t:
				lst_st2idx[s].append(0)
			lst_st2idx[s][t]=dic_st2i_val[s_val][t_val]
	return lst_st2idx

def phasedClustering(arr_exp, lst_phasedPos, metric='cosine', method='single'):
	lst_dist = scipy.spatial.distance.pdist(arr_exp, metric=metric)
	distMatrix = scipy.spatial.distance.squareform(lst_dist)
	linkageMatrix = []
	accumulatedNGene = 0
	accumulatedNGeneSubOne = 0
	totalNGene = arr_exp.shape[0]
	lst_phase = []
	dic_node2dist_nGene={}
	dic_node2phase={}
	for phase in range(len(lst_phasedPos)):
		startPos = lst_phasedPos[phase]
		if phase == len(lst_phasedPos)-1:
			endPos = arr_exp.shape[0]
		else:
			endPos = lst_phasedPos[phase+1]
		nGene = endPos - startPos
		lst_phase += [phase for i in range(nGene)]
		sub_distMatrix = distMatrix[startPos:endPos,:][:,startPos:endPos]
		sub_lst_dist = scipy.spatial.distance.squareform(sub_distMatrix)
		print phase, nGene, len(sub_lst_dist)
		sub_linkageMatrix = scipy.cluster.hierarchy.linkage(sub_lst_dist, method=method)
		maxId = 0
		maxDist = 0
		lst_subTotalID=[]
		for row in sub_linkageMatrix:
			if row[0] < nGene:
				row[0] += accumulatedNGene
			else:
				row[0] += totalNGene + accumulatedNGeneSubOne - nGene 
			if row[1] < nGene:
				row[1] += accumulatedNGene
			else:
				row[1] += totalNGene + accumulatedNGeneSubOne - nGene
			lst_subTotalID.append(row[0])
			lst_subTotalID.append(row[1])
			maxId = max(row[0],row[1],maxId)
			maxDist = max(row[2],maxDist)
			linkageMatrix.append(row)
		dic_node2dist_nGene[maxId+1]=[maxDist,nGene]
		dic_node2phase[maxId+1]=phase
		accumulatedNGene += nGene
		accumulatedNGeneSubOne += nGene - 1
	new_id = max(dic_node2dist_nGene.keys())+1
	lst_prvNode = sorted(dic_node2dist_nGene.keys(),key=lambda x:dic_node2phase[x])
	while True:
		if len(lst_prvNode) <= 1:
			break
		minI1, minI2 = None, None
		minN1, minN2 = None, None
		minDist = np.inf
		for i in range(len(lst_prvNode)-1):
			n1, n2 = lst_prvNode[i], lst_prvNode[i+1]
			dist = dic_node2dist_nGene[n1][0] + dic_node2dist_nGene[n2][0]
			if dist < minDist:
				minI1, minI2 = i, i+1
				minN1, minN2 = n1, n2
				minDist = dist
		nNode = dic_node2dist_nGene[minN1][1] + dic_node2dist_nGene[minN2][1]
		dist = dic_node2dist_nGene[minN1][0] + dic_node2dist_nGene[minN2][0]
		linkageMatrix.append([minN1,minN2,dist,nGene])
		dic_node2dist_nGene[new_id]=[dist,nNode]
		phase = min(dic_node2phase[minN1],dic_node2phase[minN2])
		dic_node2phase[new_id]=phase
		del lst_prvNode[minI2]
		del lst_prvNode[minI1]
		lst_prvNode.insert(minI1,new_id)
		new_id += 1
	return linkageMatrix

if __name__ == "__main__":
	parser=argparse.ArgumentParser(
		usage='''\
	%(prog)s [options] infile gene2phase
	example: %(prog)s infile gene2phase
	''')
	
	parser.add_argument('infile', help='DEG profile')
	parser.add_argument('gene2phase', help='gene to phase file')
	parser.add_argument('-minTimePoint', required=False, type=int, default=4, help='minimum number of time-points')
	parser.add_argument('-maxNormalization', required=False, type=bool, default=True, help='whether to perform max normalization')
	parser.add_argument('-power', required=False, type=int, default=3, help='power of emphasization (Must be an odd number to preserve the sign')
	parser.add_argument('-vmin', required=False, type=float, default=-0.4, help='the lower bound of gene expression change to be shown in the plot')
	parser.add_argument('-vmax', required=False, type=float, default=0.4, help='the upper bound of gene expression change to be shown in the plot')
	parser.add_argument('-o', required=False, metavar='str', default='result', help='outfile')
	args=parser.parse_args()
	
	# read DEG profile
	lst_sample, lst_gene, arr_exp = readData(args.infile)
	lst_st2idx = label2struct(lst_sample)
	dic_gene2idx={}
	for i, gene in enumerate(lst_gene):
		dic_gene2idx[gene]=i

	# read gene2phase
	IF=open(args.gene2phase,'r')
	dic_phase2idxs={}
	for line in IF:
		gene, phase = line.strip().split('\t')
		phase = int(phase)
		if phase not in dic_phase2idxs:
			dic_phase2idxs[phase]=[]
		dic_phase2idxs[phase].append(dic_gene2idx[gene])

	# filter out small-sized gene clusters
	set_phase=dic_phase2idxs.keys()
	for phase in set_phase:
		if len(dic_phase2idxs[phase]) <= 2:
			del dic_phase2idxs[phase]

	# prepare colors for columns (samples & time-points)
	lst_colColor = []
	lst_cidx=[]
	for s in range(len(lst_st2idx)):
		if len(lst_st2idx[s]) < args.minTimePoint:
			continue
		for t in range(len(lst_st2idx[s])):
			lst_cidx.append(lst_st2idx[s][t])
		lst_colColor+=sns.color_palette("Greys",n_colors=len(lst_st2idx[s]))

	# prepare colors for rows (genes)
	lst_rowColor = []
	lst_gidxs = []
	lst_phasePos = []
	nPhase = len(dic_phase2idxs)
	for phaseIdx, phase in enumerate(sorted(dic_phase2idxs.keys())):
		lst_phasePos.append(len(lst_gidxs))
		lst_gidxs += dic_phase2idxs[phase]
		color = sns.color_palette("rainbow_r", n_colors=nPhase, desat=0.6)[phaseIdx]
		lst_rowColor += [color for i in range(len(dic_phase2idxs[phase]))]
	lst_gene = [lst_gene[i] for i in lst_gidxs]

	# select subset of expression array
	arr_exp = arr_exp[lst_gidxs,:][:,lst_cidx]

	# max normalization
	if args.maxNormalization == True:
		start=0
		for s in range(len(lst_st2idx)):
			if len(lst_st2idx[s]) < args.minTimePoint:
				continue
			lst_cidx=[]
			for t in range(len(lst_st2idx[s])):
				lst_cidx.append(start+t)
			for i in range(len(lst_gene)):
				x = np.amax(np.abs(arr_exp[:,lst_cidx]),axis=1)
				arr_exp[i,lst_cidx] /= x[i]
			start += len(lst_st2idx[s])

	# emphasize by power
	arr_exp=np.power(arr_exp,args.power)

	# hierarchical clustering of genes with maintenance of phase
	linkageMatrix = phasedClustering(arr_exp, lst_phasePos, metric='cosine', method='single')

	# make heatmap plot
	cm=sns.clustermap(arr_exp, cmap='RdYlBu_r', row_cluster=True, col_cluster=False, row_linkage=linkageMatrix, col_colors=lst_colColor, row_colors=lst_rowColor, xticklabels=False, yticklabels=False, vmin=args.vmin, vmax=args.vmax)

	# show plot
	plt.show()
