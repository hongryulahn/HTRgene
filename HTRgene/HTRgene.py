#!/usr/bin/env python
import sys
import argparse
import numpy as np
import os
import itertools
import time
import random
from multiprocessing import Pool, Lock
import scipy.stats
import scipy.spatial.distance
import subprocess
import statsmodels.sandbox.stats.multicomp

from matplotlib import pyplot as plt

def ranking(x, axis=0, method='average', reverse=False):
	"""Ranking numpy array (1D or 2D dimension).

	Return ranking of numpy array.
	If x = [[1,2,3],
		[4,3,2],
		[6,5,1]]
	it will be returned 	[[0,0,2],
				 [1,1,1],
				 [2,2,0]]
	Args:
		x: numpy array
		axis: Ranking for each column (0) or for each row (1) when x is 2D-dimensional array
			(default = 0)
		method: How to assign ranking for the same value (see scipy.stats.rankdata)
			(default = 'average')
		reverse: Decending order
			(default = False)
	Returns:
		A ranking 2D numpy array
	Raises:
		pass
	"""
	
	b=np.array(x)
	if reverse:
		b=-1.0*b
	nDimension = len(b.shape)
	if nDimension == 1:
		return scipy.stats.rankdata(b, method=method)
	elif nDimension == 2:
		nRow, nCol = b.shape
		if axis == 0:
			for col in range(nCol):
				b[:,col] = scipy.stats.rankdata(b[:,col],method=method)
		elif axis == 1:
			for row in range(nRow):
				b[row,:] = scipy.stats.rankdata(b[row,:],method=method)
		return b
	else:
		print('I don\' know how to rank %d dimensional data'%(nDimension))
		sys.exit(1)

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
	dic_str2i_val={}
	for i, sample in enumerate(lst_sample):
		s_val, t_val, r_val = sample.split('_')
		s_val, t_val, r_val = int(s_val.strip('s')), int(t_val.strip('t')), int(r_val.strip('r'))
		if s_val not in dic_str2i_val:
			dic_str2i_val[s_val]={}
		if t_val not in dic_str2i_val[s_val]:
			dic_str2i_val[s_val][t_val]={}
		if r_val not in dic_str2i_val[s_val][t_val]:
			dic_str2i_val[s_val][t_val][r_val]=i

	lst_str2idx=[]
	for s, s_val in enumerate(sorted(dic_str2i_val.keys(),key=lambda x:float(x))):
		for t, t_val in enumerate(sorted(dic_str2i_val[s_val].keys(),key=lambda x:float(x))):
			for r, r_val in enumerate(sorted(dic_str2i_val[s_val][t_val].keys(),key=lambda x:float(x))):
				if len(lst_str2idx) <= s:
					lst_str2idx.append([])
				if len(lst_str2idx[s]) <= t:
					lst_str2idx[s].append([])
				if len(lst_str2idx[s][t]) <= r:
					lst_str2idx[s][t].append(0)
				lst_str2idx[s][t][r]=dic_str2i_val[s_val][t_val][r_val]
	return lst_str2idx

'''
def arrToSarr(arr_exp, lst_sample):
	dic_etr2column={}
	for i, sample in enumerate(lst_sample):
		e, t, r = sample.split('_')
		e, t, r = int(e.strip('e')), int(t.strip('t')), int(r.strip('r'))
		if e not in dic_etr2column:
			dic_etr2column[e]={}
		if t not in dic_etr2column[e]:
			dic_etr2column[e][t]={}
		if r not in dic_etr2column[e][t]:
			dic_etr2column[e][t][r]=i

	print arr_exp.shape
	sarr_exp_gstr=[None for g in range(arr_exp.shape[0])]
	for g in range(arr_exp.shape[0]):
		sarr_exp_gstr[g]=[None for e in range(len(dic_etr2column))]
		for e, e_val in enumerate(sorted(dic_etr2column.keys(),key=lambda x:float(x))):
			sarr_exp_gstr[g][e]=[None for t in range(len(dic_etr2column[e_val]))]
			for t, t_val in enumerate(sorted(dic_etr2column[e_val].keys(),key=lambda x:float(x))):
				sarr_exp_gstr[g][e][t]=[None for r in range(len(dic_etr2column[e_val][t_val]))]
				for r, r_val in enumerate(sorted(dic_etr2column[e_val][t_val].keys(),key=lambda x:float(x))):
					#print g, e, t, r, dic_etr2column[e_val][t_val][r_val]
					#print arr_exp[g,:]
					sarr_exp_gstr[g][e][t][r]=arr_exp[g,dic_etr2column[e_val][t_val][r_val]]
		print g
	return sarr_exp_gstr
'''
	
def normailze(arr_exp, lst_str2idx, platform='microarray'):
	arr_basement = np.zeros((arr_exp.shape[0],len(lst_str2idx)))
	for s in range(len(lst_str2idx)):
		lst_idx_t0=lst_str2idx[s][0]
		if platform=='microarray':
			arr_basement[:,s]=np.mean(arr_exp[:,lst_idx_t0],axis=1)
		elif platform=='RNA-seq':
			arr_basement[:,s]=np.mean(np.log(arr_exp[:,lst_idx_t0]+1),axis=1)
		
	y=np.zeros(arr_exp.shape)
	for s in range(len(lst_str2idx)):
		for t in range(len(lst_str2idx[s])):
			for r in range(len(lst_str2idx[s][t])):
				i = lst_str2idx[s][t][r]
				if platform=='microarray':
					y[:,i] = arr_exp[:,i] - arr_basement[:,s]
				elif platform=='RNA-seq':
					y[:,i] = np.log(arr_exp[:,i]+1) - arr_basement[:,s]
	return y

def profileDEG(arr_exp, lst_str2idx, lst_gene=None, valueType='absolute_binary',passOut=None,outdir=None):
	if outdir == None:
		tmp1='tmp1'
		tmp2='tmp2'
		tmp3='tmp3'
		tmp4='tmp4'
	else:
		tmp1=outdir+'/tmp1'
		tmp2=outdir+'/tmp2'
		tmp3=outdir+'/tmp3'
		tmp4=outdir+'/tmp4'
	if passOut == None:
		OF1=open(tmp1,'w')
		OF2=open(tmp2,'w')
		OF3=open(tmp3,'w')
		lst_label=[]
		for s in range(len(lst_str2idx)):
			for t in range(len(lst_str2idx[s])):
				label_st = '%d_%d'%(s,t)
				for r in range(len(lst_str2idx[s][t])):
					label_str = '%d_%d_%d'%(s,t,r)
					OF2.write(label_str+'\t'+label_st+'\n')
					lst_label.append(label_str)
				if t != 0:
					label_s0 = '%d_%d'%(s,0)
					OF3.write(label_s0+'\t'+label_st+'\n')
		OF1.write('ID'+'\t'+'\t'.join(lst_label)+'\n')
		
		for g in range(arr_exp.shape[0]):
			if lst_gene == None:
				name = 'g%09d'%(g)
			else:
				name = lst_gene[g]
			OF1.write(name+'\t'+'\t'.join(map(str,arr_exp[g,:]))+'\n')
		OF1.close()
		OF2.close()
		OF3.close()
		subprocess.call(["Rscript", "DEGlimma.R", tmp1, tmp2, tmp3, valueType, tmp4])
	IF = open(tmp4,'r')
	lst_st2idx=[]
	lst_label=IF.readline().strip().split('\t')[1:]
	for i, label in enumerate(lst_label):
		l0,l1 = label.split(',')
		s=int(l1.split('_')[0])
		t=int(l1.split('_')[1])-1
		if len(lst_st2idx) <= s:
			lst_st2idx.append([])
		if len(lst_st2idx[s]) <= t:
			lst_st2idx[s].append(0)
		lst_st2idx[s][t]=i
	lst_DEGs=[]
	for line in IF:
		s=line.strip().split('\t')
		name, DEGs = s[0], map(float,s[1:])
		lst_DEGs.append(DEGs)
	y=np.array(lst_DEGs)
	IF.close()
	return y, lst_st2idx

def mergeDEGprofile(arr_DEG, lst_st2idx):
	y=np.zeros((arr_DEG.shape[0], len(lst_st2idx)))
	for s in range(len(lst_st2idx)):
		lst_i = lst_st2idx[s]
		y[:,s] = np.max(arr_DEG[:,lst_i],axis=1)
	return y

def testDEGprofile(arr_DEG, method='sampling', nTrial = 1000000, seed=None):
	if method == 'sampling':
		if seed != None:
			np.random.seed(seed)
		nGene = arr_DEG.shape[0]
		#nTrial = nGene * 100
		lst_randomVal = []
		x = arr_DEG.copy()
		while True:
			if len(lst_randomVal) >= nTrial:
				break
			for col in range(x.shape[1]):
				x[:,col]=np.random.permutation(x[:,col])
			for row in range(x.shape[0]):
				if len(lst_randomVal) >= nTrial:
					break
				lst_randomVal.append(np.sum(x[row,:]))
			#lst_randomVal.append(np.sum([x[row,col] for col, row in enumerate(np.random.choice(nGene, arr_DEG.shape[1]))]))
		lst_pval = np.ones(nGene)
		lst_gene_val = sorted([(g, np.sum(arr_DEG[g,:])) for g in range(nGene)],key=lambda x:x[1], reverse=True)
		lst_randomVal.sort(reverse=True)
		lst_randomVal.append(-np.inf)
		i,j=0,0
		while True:
			if i >= nGene:
				break
			if lst_randomVal[j] == lst_randomVal[j+1]:
				j+=1
				continue
			gene, val = lst_gene_val[i]
			randomVal = lst_randomVal[j]
			if val < randomVal:
				j+=1
			elif val >= randomVal:
				lst_pval[gene] = float(j)/nTrial
				i+=1
	return lst_pval

def testCluster(distMatrix, dic_cluster2idxs, method='sampling', nTrial = 10000):
	if method == 'sampling':
		nGene = distMatrix.shape[0]
		lst_cluster = sorted(dic_cluster2idxs.keys(),key=lambda x: len(dic_cluster2idxs[x]))
		max_size = len(dic_cluster2idxs[lst_cluster[-1]])

def selectSignificant(lst_pval, method='pcut', multipleTestCorrection='bonferroni', alpha=0.05, topN=None):
	lst_idx=[]
	if method == 'top':
		lst_ranking = ranking(lst_pval)
		for i in range(len(lst_pval)):
			if lst_ranking[i] <= topN:
				lst_idx.append(i)
	if method == 'pcut':
		if multipleTestCorrection != None:
			reject,lst_pval,alphacSidak,alphacBonf = statsmodels.sandbox.stats.multicomp.multipletests(lst_pval,method=multipleTestCorrection,alpha=alpha)
		for i in range(len(lst_pval)):
			if lst_pval[i] <= alpha:
				lst_idx.append(i)
	return lst_idx

def getLowerBoundDistance(distMatrix, alpha=0.05):
	sortedDistMatrix=np.sort(distMatrix, axis=1)
	nGene = sortedDistMatrix.shape[0]
	lst_lowerBoundDistance=[]
	lst_lowerBoundDistance.append(None)
	lst_lowerBoundDistance.append(None)
	for size in range(2,nGene):
		lst_dist = np.sort(sortedDistMatrix[:,size-1])
		lst_lowerBoundDistance.append(lst_dist[int(nGene*alpha)])
	return lst_lowerBoundDistance

def skmeansWrapper(arr_exp, k, passOut=None, seed=None, outdir=None):
	if outdir == None:
		tmp5='tmp5'
		tmp6='tmp6'
	else:
		tmp5=outdir+'/tmp5'
		tmp6=outdir+'/tmp6'
	if passOut == None:
		OF5=open(tmp5,'w')
		lst_label=[str(i) for i in range(arr_exp.shape[1])]
		OF5.write('ID'+'\t'+'\t'.join(lst_label)+'\n')
		for g in range(arr_exp.shape[0]):
			name = 'g%09d'%(g)
			OF5.write(name+'\t'+'\t'.join(map(str,arr_exp[g,:]))+'\n')
		OF5.close()
		if seed == None:
			subprocess.call(["Rscript", "run_skmeans.R", tmp5, str(k), 'None', tmp6])
		else:
			subprocess.call(["Rscript", "run_skmeans.R", tmp5, str(k), str(seed), tmp6])
	IF = open(tmp6,'r')
	dic_cluster2idxs={}
	for line in IF:
		s=line.strip().split('\t')
		name, cluster = s[0], int(s[1])
		idx=int(name.replace('g',''))
		if cluster not in dic_cluster2idxs:
			dic_cluster2idxs[cluster]=set()
		dic_cluster2idxs[cluster].add(idx)
	IF.close()
	return dic_cluster2idxs

def testCluster(dic_cluster2idxs, distMatrix, method='sampling', nTrial=1000, seed=None):
	if seed != None:
		np.random.seed(seed)
	lst_cid2cluster=[]
	dic_cluster2cid={}
	lst_cid_position=[]
	lst_idxs=[]
	for cid, (cluster, idxs) in enumerate(sorted(dic_cluster2idxs.items())):
		lst_cid2cluster.append(cluster)
		dic_cluster2cid[cluster]=cid
		lst_cid_position.append(len(lst_idxs))
		lst_idxs += sorted(idxs)
	lst_cid_position.append(len(lst_idxs))
	original_idxs=np.array(lst_idxs)
	lst_meanDist=[]
	lst_dist=[]
	for cid in range(nCluster):
		tmpIdxs = original_idxs[lst_cid_position[cid]:lst_cid_position[cid+1]]
		if len(tmpIdxs) <= 1:
			lst_meanDist.append(0.0)
			lst_dist.append(0.0)
			continue
		tmpDistMatrix = distMatrix[tmpIdxs,:][:,tmpIdxs]
		center = sorted([(i,np.sum(tmpDistMatrix[i,:])) for i in range(len(tmpIdxs))],key=lambda x:x[1])[0][0]
		lst_meanDist.append(np.mean(tmpDistMatrix[center,:]))
		lst_dist+=list(tmpDistMatrix[center,:])
	nCluster=len(dic_cluster2idxs)
	nIdx=len(lst_idxs)
	lst_randomMeanDist=[[] for i in range(nCluster)]
	lst_randomDist=[[] for i in range(nCluster)]
	for t in nTrial:
		random_idxs=np.random.permutation(nIdx)
		for cid in range(nCluster):
			tmpIdxs = random_idxs[lst_cid_position[cid]:lst_cid_position[cid+1]]
			if len(tmpIdxs) <= 1:
				lst_meanDist[cid].append(0.0)
				lst_dist[cid].append(0.0)
				continue
			tmpDistMatrix = distMatrix[tmpIdxs,:][:,tmpIdxs]
			center = sorted([(i,np.sum(tmpDistMatrix[i,:])) for i in range(len(tmpIdxs))],key=lambda x:x[1])[0][0]
			lst_randomMeanDist[cid].append(np.mean(tmpDistMatrix[center,:]))
			lst_randomDist[cid]+=list(tmpDistMatrix[center,:])
	#for cid in range(nCluster):
	return lst_clusterPval, lst_idxPval 

def getResponseTimeVector(arr_exp, lst_str2idx):
	lst_posBreakpoint=[]
	lst_posAllTval=[]
	lst_negBreakpoint=[]
	lst_negAllTval=[]
	for s in range(len(lst_str2idx)):
		lst_timepointTval=[]
		for tidx in range(1,len(lst_str2idx[s])):
			beforeIdx=[]
			afterIdx=[]
			for tidx2 in range(len(lst_str2idx[s])):
				for r, idx in enumerate(lst_str2idx[s][tidx2]):
					if tidx2 < tidx:
						beforeIdx.append(idx)
					else:
						afterIdx.append(idx)
			if len(beforeIdx) <= 1 or len(afterIdx) <= 1:
				#tval, pval = 0.0, 1.0
				mean_tval = 0.0
			else:
				beforeExp=arr_exp[:,beforeIdx]
				afterExp=arr_exp[:,afterIdx]
				lst_tval=[]
				for gidx in range(arr_exp.shape[0]):
					tval, pval = scipy.stats.ttest_ind(afterExp[gidx,:],beforeExp[gidx,:])
					lst_tval.append(tval)
				mean_tval = np.mean(lst_tval)
			lst_timepointTval.append([tidx,mean_tval])
		breakpoint, tval_breakpoint = sorted(lst_timepointTval,key=lambda x:x[1],reverse=True)[0]
		lst_posBreakpoint.append(breakpoint)
		lst_posAllTval.append(tval_breakpoint)
		breakpoint, tval_breakpoint = sorted(lst_timepointTval,key=lambda x:x[1])[0]
		lst_negBreakpoint.append(breakpoint)
		lst_negAllTval.append(tval_breakpoint)
	if abs(np.mean(lst_posAllTval)) > abs(np.mean(lst_negAllTval)):
		lst_breakpoint=lst_posBreakpoint
		allTval=np.mean(lst_posAllTval)
	else:
		lst_breakpoint=lst_negBreakpoint
		allTval=np.mean(lst_negAllTval)
	return lst_breakpoint, allTval

def getResponseTimeVector2(arr_exp, lst_str2idx):
	lst_posBreakpoint=[]
	lst_posTval=[]
	lst_negBreakpoint=[]
	lst_negTval=[]
	for s in range(len(lst_str2idx)):
		lst_timepointTval=[]
		for tidx in range(1,len(lst_str2idx[s])):
			beforeIdx=[]
			afterIdx=[]
			for tidx2 in range(len(lst_str2idx[s])):
				for r, idx in enumerate(lst_str2idx[s][tidx2]):
					if tidx2 < tidx:
						beforeIdx.append(idx)
					else:
						afterIdx.append(idx)
			if len(beforeIdx) <= 1 or len(afterIdx) <= 1:
				#tval, pval = 0.0, 1.0
				mean_tval = 0.0
			else:
				beforeExp=arr_exp[:,beforeIdx].flatten()
				afterExp=arr_exp[:,afterIdx].flatten()
				tval, pval = scipy.stats.ttest_ind(afterExp,beforeExp)
			lst_timepointTval.append([tidx,tval])
		breakpoint, tval = sorted(lst_timepointTval,key=lambda x:x[1],reverse=True)[0]
		lst_posBreakpoint.append(breakpoint)
		lst_posTval.append(tval)
		breakpoint, tval = sorted(lst_timepointTval,key=lambda x:x[1])[0]
		lst_negBreakpoint.append(breakpoint)
		lst_negTval.append(tval)
	if abs(np.mean(lst_posTval)) >= abs(np.mean(lst_negTval)):
		flag_up = True
		lst_breakpoint=lst_posBreakpoint
	else:
		flag_up = False
		lst_breakpoint=lst_negBreakpoint
	lst_final_breakpoint, lst_final_tval = [], []
	for s in range(len(lst_str2idx)):
		breakpoint = lst_breakpoint[s]
		beforeIdx=[]
		afterIdx=[]
		for tidx2 in range(len(lst_str2idx[s])):
			for r, idx in enumerate(lst_str2idx[s][tidx2]):
				if tidx2 < breakpoint:
					beforeIdx.append(idx)
				else:
					afterIdx.append(idx)
		if len(beforeIdx) <= 1 or len(afterIdx) <= 1:
			#tval, pval = 0.0, 1.0
			mean_tval = 0.0
		else:
			beforeExp=arr_exp[:,beforeIdx]
			afterExp=arr_exp[:,afterIdx]
			lst_tval=[]
			for gidx in range(arr_exp.shape[0]):
				tval, pval = scipy.stats.ttest_ind(afterExp[gidx,:],beforeExp[gidx,:])
				lst_tval.append(tval)
			mean_tval = np.mean(lst_tval)
		leftBreakpoint, leftMeanTval = breakpoint, mean_tval
		while True:
			if leftBreakpoint <= 1:
				break
			leftBreakpoint -= 1
			beforeIdx=[]
			afterIdx=[]
			for tidx2 in range(len(lst_str2idx[s])):
				for r, idx in enumerate(lst_str2idx[s][tidx2]):
					if tidx2 < leftBreakpoint:
						beforeIdx.append(idx)
					else:
						afterIdx.append(idx)
			if len(beforeIdx) <= 1 or len(afterIdx) <= 1:
				#tval, pval = 0.0, 1.0
				mean_tval = 0.0
			else:
				beforeExp=arr_exp[:,beforeIdx]
				afterExp=arr_exp[:,afterIdx]
				lst_tval=[]
				for gidx in range(arr_exp.shape[0]):
					tval, pval = scipy.stats.ttest_ind(afterExp[gidx,:],beforeExp[gidx,:])
					lst_tval.append(tval)
				mean_tval = np.mean(lst_tval)
			if flag_up:
				if mean_tval > leftMeanTval:
					leftMeanTval = mean_tval
				else:
					leftBreakpoint += 1
					break
			else:
				if mean_tval < leftMeanTval:
					leftMeanTval = mean_tval
				else:
					leftBreakpoint += 1
					break
		rightBreakpoint, rightMeanTval = breakpoint, mean_tval
		while True:
			if rightBreakpoint >= len(lst_str2idx[s])-1:
				break
			rightBreakpoint += 1
			beforeIdx=[]
			afterIdx=[]
			for tidx2 in range(len(lst_str2idx[s])):
				for r, idx in enumerate(lst_str2idx[s][tidx2]):
					if tidx2 < rightBreakpoint:
						beforeIdx.append(idx)
					else:
						afterIdx.append(idx)
			if len(beforeIdx) <= 1 or len(afterIdx) <= 1:
				#tval, pval = 0.0, 1.0
				mean_tval = 0.0
			else:
				beforeExp=arr_exp[:,beforeIdx]
				afterExp=arr_exp[:,afterIdx]
				lst_tval=[]
				for gidx in range(arr_exp.shape[0]):
					tval, pval = scipy.stats.ttest_ind(afterExp[gidx,:],beforeExp[gidx,:])
					lst_tval.append(tval)
				mean_tval = np.mean(lst_tval)
			if flag_up:
				if mean_tval > rightMeanTval:
					rightMeanTval = mean_tval
				else:
					rightBreakpoint -= 1
					break
			else:
				if mean_tval < rightMeanTval:
					rightMeanTval = mean_tval
				else:
					rightBreakpoint -= 1
					break
		if flag_up > 0:
			if leftMeanTval > rightMeanTval:
				breakpoint = leftBreakpoint
				meanTval = leftMeanTval
			else:
				breakpoint = rightBreakpoint
				meanTval = rightMeanTval
		else:
			if leftMeanTval < rightMeanTval:
				breakpoint = leftBreakpoint
				meanTval = leftMeanTval
			else:
				breakpoint = rightBreakpoint
				meanTval = rightMeanTval
		lst_final_breakpoint.append(breakpoint)
		lst_final_tval.append(meanTval)
	return lst_final_breakpoint, np.mean(lst_final_tval)

def breakpointOrder(breakpoint1, breakpoint2):
	lessOrEqual=[]
	greaterOrEqual=[]
	for b1,b2 in zip(breakpoint1, breakpoint2):
		if b1 == None or b2 == None:
			lessOrEqual.append(True)
			greaterOrEqual.append(True)
		elif b1 < b2:
			lessOrEqual.append(True)
			greaterOrEqual.append(False)
		elif b1 > b2:
			lessOrEqual.append(False)
			greaterOrEqual.append(True)
		elif b1 == b2:
			lessOrEqual.append(True)
			greaterOrEqual.append(True)
	if sum(lessOrEqual) == len(breakpoint1) and sum(greaterOrEqual) == len(breakpoint1):
		return '=='
	elif sum(lessOrEqual) == len(breakpoint1):
		return '<'
	elif sum(greaterOrEqual) == len(breakpoint1):
		return '>'
	else:
		return 'conflict'

def orderingCluster(dic_cluster2rt_tstats):
	lst_ordered_breakpoint=[]
	dic_ordered_breakpoint2clusters={}
	for cluster in sorted(dic_cluster2rt_tstats.keys(), key=lambda x:abs(dic_cluster2rt_tstats[x][1]), reverse=True):
		breakpoint, tstats = dic_cluster2rt_tstats[cluster]
		breakpoint = tuple(breakpoint)
		if len(lst_ordered_breakpoint) == 0:
			lst_ordered_breakpoint.append(breakpoint)
			dic_ordered_breakpoint2clusters[breakpoint]=set([cluster])
			continue
		if breakpoint in dic_ordered_breakpoint2clusters:
			dic_ordered_breakpoint2clusters[breakpoint].add(cluster)
			continue
		for i, cur_breakpoint in enumerate(lst_ordered_breakpoint):
			order = breakpointOrder(breakpoint,cur_breakpoint)
			if order == '<':
				lst_ordered_breakpoint.insert(i,breakpoint)
				dic_ordered_breakpoint2clusters[breakpoint]=set([cluster])
				break
			elif order == '==':
				if breakpoint == cur_breakpoint:
					continue
				dic_ordered_breakpoint2clusters[breakpoint].add(cluster)
				break
			elif order == '>':
				if i == len(lst_ordered_breakpoint)-1:
					lst_ordered_breakpoint.append(breakpoint)
					dic_ordered_breakpoint2clusters[breakpoint]=set([cluster])
			elif order == 'conflict':
				break
	
	dic_cluster2order={}
	for i, breakpoint in enumerate(lst_ordered_breakpoint,1):
		set_cluster=sorted(dic_ordered_breakpoint2clusters[breakpoint])
		for cluster in set_cluster:
			dic_cluster2order[cluster]=i
	return dic_cluster2order


def run(lst_sample, lst_gene, arr_exp, DEGpcut, cosinePcut, minSizePcut, step2_k='auto', platform='microarray', seed=None, outdir=None):
	dic_gene2phase={}
	########################
	#        STEP 1        #
	########################
	# Convert array to sarr
	lst_str2idx = label2struct(lst_sample)
	#sarr_exp_gstr = arrToSarr(arr_exp, lst_sample)
	'''
		A sarr_exp (structured array of expression) has four dimensions.
		The first dimension elements, sarr_exp[.], are genes.
		The second dimension elements, sarr_exp[.][.], are samples.
		The third dimension elements, sarr_exp[.][.][.], are time-points.
		The fourth dimension elements, sarr_exp[.][.][.][.], are replicates.
		Thus, the sarr_exp[g][s][t][r] value is the expression level of gene g in sample s, timepoint t, replicate r.
	'''
	# Normalize expression
	arr_exp = normailze(arr_exp, lst_str2idx, platform=platform)

	# Profile DEG for t0 vs. tx for each sample
	arr_DEG, lst_st2idx = profileDEG(arr_exp, lst_str2idx, lst_gene=lst_gene, valueType='absolute_binary', passOut=None, outdir=outdir)

	# Merge DEG profile for each sample
	arr_DEG = mergeDEGprofile(arr_DEG, lst_st2idx)
	'''
		A arr_DEGprofile is a numpy array of two dimensions.
		That is, the value of arr_DEGprofile[g,e] is 1(DEG) or 0(non-DEG) for gene g in sample s.
	'''
	# Compute adjusted p-value of DEG frequencies against the random frequency distribution.
	lst_frequency = np.sum(arr_DEG, axis=1)
	lst_pval = testDEGprofile(arr_DEG, method='sampling', nTrial=200000, seed=seed)
	#lst_DEGidx = selectSignificant(lst_pval, method='pcut', multipleTestCorrection='bonferroni', alpha=0.05)
	lst_DEGidx = selectSignificant(lst_pval, method='pcut', multipleTestCorrection='fdr_bh', alpha=DEGpcut)

	for gene in set(lst_gene)-set([lst_gene[i] for i in lst_DEGidx]):
		dic_gene2phase[gene]='outByDEGfreq'
	#for i in range(len(lst_gene)):
	#	print lst_gene[i], lst_pval[i]
	lst_sortedDEGidx = sorted(lst_DEGidx, key=lambda x:lst_frequency[x], reverse=True)
	lst_DEGgene = [lst_gene[i] for i in lst_sortedDEGidx]
	arr_DEGexp = arr_exp[lst_sortedDEGidx,:]
	set_pseudoRefGene = set(lst_DEGgene[0:int(0.1*len(lst_DEGgene))])

	if step2_k == 'auto':
		maxK = min(500,len(arr_DEGexp))
		lst_k = list(range(50,maxK,50))
	else:
		lst_k = [int(step2_k)]
	distMatrix=scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(arr_DEGexp, metric='cosine'))
	#lst_lowerBoundDistance = getLowerBoundDistance(distMatrix, cosinePcut)
	dic_k2gene_phase={}
	dic_k2f1={}
	for k in lst_k:
		########################
		#        STEP 2        #
		########################
		dic_cluster2idxs = skmeansWrapper(arr_DEGexp, k, seed=seed, outdir=outdir) 
		lst_clusterPval, lst_idxPval = testCluster(dic_cluster2idxs, distMatrix, method='sampling', nTrial=1000)
		'''
		set_outlier = set()
		for cluster, set_idx in dic_cluster2idxs.items():
			center = sorted([(i,np.sum(distMatrix[i,sorted(set_idx)])) for i in set_idx],key=lambda x:x[1])[0][0]
			for idx in set_idx:
				if distMatrix[center,idx] > lst_lowerBoundDistance[len(set_idx)]:
					set_outlier.add(idx)
			dic_cluster2idxs[cluster] = set([idx for idx in set_idx if idx not in set_outlier])

		for i in set_outlier:
			gene = lst_DEGgene[i]
			dic_gene2phase[gene]='outByCosDist'

		set_remove = set()
		lst_nGene = [len(set_idx) for set_idx in dic_cluster2idxs.values()]
		lst_pval=[]
		for cluster, set_idx in sorted(dic_cluster2idxs.items()):
			sizePval = scipy.stats.binom.cdf(len(set_idx),len(arr_DEGexp),1.0/float(k))
			lst_pval.append(sizePval)
			#if tstats < 0 and sizePval < minSizePcut:
			#	set_remove.add(cluster)

		reject,lst_qval,alphacSidak,alphacBonf = statsmodels.sandbox.stats.multicomp.multipletests(lst_pval,method='fdr_bh',alpha=0.05)
		#lst_qval = selectSignificant(lst_pval, method='pcut', multipleTestCorrection='fdr_bh', alpha=DEGpcut)
		for i, (cluster, set_idx) in enumerate(sorted(dic_cluster2idxs.items())):
			tstats=lst_tstats[i]
			pval = lst_pval[i]
			sizePval = lst_qval[i]
			print cluster, len(set_idx), tstats, pval, sizePval, np.mean(lst_nGene), len(arr_DEGexp)*1.0/float(k)
			if tstats > 0 and sizePval < minSizePcut:
				set_remove.add(cluster)

		for cluster in set_remove:
			for idx in dic_cluster2idxs[cluster]:
				gene=lst_DEGgene[idx]
				dic_gene2phase[gene]='outBySmallCluster'
			del dic_cluster2idxs[cluster]
		'''
		
		########################
		#        STEP 3        #
		########################
		dic_cluster2rt_tstats={}
		for i, (cluster, gidxs) in enumerate(dic_cluster2idxs.items()):
			#if i > 3: break
			rt, tstats = getResponseTimeVector(arr_DEGexp[sorted(gidxs),:], lst_str2idx)
			dic_cluster2rt_tstats[cluster]=[rt,tstats]
			#print cluster, rt, tstats

		########################
		#        STEP 4        #
		########################
		dic_gene2phase={}
		dic_cluster2phase = orderingCluster(dic_cluster2rt_tstats)
		for cluster, phase in dic_cluster2phase.items():
			for i in dic_cluster2idxs[cluster]:
				gene = lst_DEGgene[i]
				dic_gene2phase[gene]=phase

		dic_k2gene_phase[k]=dic_gene2phase
		set_gene=set(dic_gene2phase.keys())
		TP = len(set_pseudoRefGene & set_gene)
		FP = len(set_gene - set_pseudoRefGene)
		FN = len(set_pseudoRefGene - set_gene)
		TN = len(lst_DEGgene) - TP - FP - FN
		if TP == 0:
			f1 = 0.0
		else:
			precision = float(TP)/(TP+FP)
			recall = float(TP)/(TP+FN)
			f1 = 2*(precision*recall)/(precision+recall)
		dic_k2f1[k]=f1
		print k, TP, FP, FN, TN, f1
	k, f1 = sorted(dic_k2f1.items(),key=lambda x:x[1], reverse=True)[0]
	return dic_k2gene_phase[k]



if __name__ == "__main__":
	parser=argparse.ArgumentParser(
		usage='''\
	%(prog)s [options] expfile -o out.txt
	example: %(prog)s expfile -o out.txt
	''')
	
	parser.add_argument('expfile', help='Gene expression file')
	parser.add_argument('-step2_k', required=False, default='auto', help='The number of clusters for Step 2')
	parser.add_argument('-DEGpcut', required=False, type=float, default=0.05, help='')
	parser.add_argument('-cosinePcut', required=False, type=float, default=0.05, help='DEGratio cutoff')
	parser.add_argument('-minSizePcut', required=False, type=float, default=0.01, help='ttest|MI')
	parser.add_argument('-seed', required=False, type=int, default=20170211, help='Seed value to fix random process')
	parser.add_argument('-o', required=False, metavar='str', default='result', help='out directory')
	args=parser.parse_args()
	
	if not os.path.exists(args.o):
		os.makedirs(args.o)

	lst_sample, lst_gene, arr_exp = readData(args.expfile)

	#dic_orderedCluster2genes = run(lst_sample=lst_sample, lst_gene=lst_gene, arr_exp=arr_exp, DEGpcut=args.DEGpcut, cosinePcut=args.cosinePcut, minSizePcut=args.minSizePcut, step2_k='auto', platform='microarray')
	dic_gene2phase = run(lst_sample=lst_sample, lst_gene=lst_gene, arr_exp=arr_exp, DEGpcut=args.DEGpcut, cosinePcut=args.cosinePcut, minSizePcut=args.minSizePcut, step2_k=args.step2_k, platform='microarray', seed=args.seed, outdir=args.o)

	OF=open(args.o+'/gene2phase.txt','w')
	for gene, phase in dic_gene2phase.items():
		OF.write('\t'.join(map(str,[gene, phase]))+'\n')
	#for oc, lst_gene in dic_orderedCluster2genes.items():
	#	for gene in lst_gene:
	#		OF.write('\t'.join(map(str,[gene, oc]))+'\n')
