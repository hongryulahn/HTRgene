"""
HTRgene

Created by Hongryul Ahn on 2018-04-07.
Copyright (c) 2018 Hongryul Ahn. All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

import sys
import argparse
import numpy as np
import os
import itertools
import random
import scipy.stats
import subprocess
import statsmodels.sandbox.stats.multicomp

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

def isNumber(x):
	try:
		float(x)
		return True
	except:
		return False

def readData(expfiles):
	"""Read data from expfile.

	Return a list of samples, a list of genes, and gene expression array (numpy.array).

	Args:
		expfiles: gene expression matrix files
			An expfile includes multiple time-series samples
			Each sample ID in header line of the file must follow the format SID_TID_RID
				SID is a string to indicates the time-series sample
				TID is a number to indicates the time point (each sample must have at least one 0 time point)
				RID is a number to indicates the biological replicate
			The SID, TID, RID are sorted and assigned a new interger ID starting from 0.
	Returns:
		a list of samples, a list of genes, and gene expression array
	Raises:
		pass
	"""

	dic_gene2exps={}
	lst_newSample = []
	countSID=0
	for expfile in expfiles:
		IF=open(expfile,'r')
		line=IF.readline()
		lst_sample=line.strip().split('\t')[1:]

		lst_tag=[]
		dic_SID2newSID={}
		for idx, sample in enumerate(lst_sample):
			SID, TID, RID = sample.split('_')
			assert(isNumber(TID)),"TID must be a number"
			assert(isNumber(RID)),"RID must be a number"
			if SID not in dic_SID2newSID:
				dic_SID2newSID[SID]=countSID
				countSID+=1
			newSID = dic_SID2newSID[SID]
			lst_tag.append([sample, newSID, float(TID), float(RID), idx])
		lst_tag.sort(key=lambda x:(x[1],x[2],x[3]))
		lst_idxMap = []
		prvSID = None
		for sample, SID, TID, RID, idx in lst_tag:
			if SID != prvSID:
				assert(TID == 0.0),"A sample must have at least one 0 time point"
				prvTID = None
				newTID = -1
			if TID != prvTID:
				newTID += 1
				newRID = 0
			lst_newSample.append('%d_%d_%d'%(SID,newTID,newRID))
			lst_idxMap.append(idx)
			prvSID = SID
			prvTID = TID
			newRID += 1
		lst_exps=[]
		for line in IF:
			s=line.strip().split('\t')
			gene, lst_exp = s[0], [float(s[idx+1]) for idx in lst_idxMap]
			if gene not in dic_gene2exps:
				dic_gene2exps[gene] = []
			dic_gene2exps[gene]+=lst_exp

	lst_gene = []
	lst_exps = []
	for gene, lst_exp in sorted(dic_gene2exps.items()):
		if len(lst_exp) == len(lst_newSample):
			lst_gene.append(gene)
			lst_exps.append(lst_exp)
	arr_exp = np.array(lst_exps)
	return lst_newSample, lst_gene, arr_exp

def label2struct(lst_sample):
	"""Bulid a sample-timepoint-replicate structure from sample labels.

	Args:
		lst_sample: a list of labels of samples whose format is SID_TID_RID.
			For instance, 1_0_1 is a possible label where SID(sample ID), TID(timepoint), RID(replicate ID) = 1, 0, 1.
			SID, TID and RID must be integers.
			The first index of SID, TID, RID must start by 0.
			The TID = 0 (the zero timepoint) means non-stressed sample.
	Returns:
		lst_str2idx: a list of mapping sample, time-point, and replicate to column index.
			lst_str2idx[s] returns the list of time-points of sample s.
			lst_str2idx[s][t] returns the list of replicates of time-points t of sample s.
			lst_str2idx[s][t][r] returns the column index of replicate r of time-point t of sample s.
	Raises:
		pass
	"""
	dic_str2i_val={}
	for i, sample in enumerate(lst_sample):
		s_val, t_val, r_val = sample.split('_')
		t_val = int(t_val)
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

def normailze(arr_exp, lst_str2idx, platform='microarray'):
	"""Bulid a sample-timepoint-replicate structure from sample labels.

	Args:
		expfile: gene expression matrix file
	Returns:

		a list of samples, a list of genes, and gene expression array
	Raises:
		pass
	"""

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

def profileDEG(arr_exp, lst_str2idx, lst_gene=None, valueType='absolute_binary',passBy=False,outdir=None):
	"""profiling DEGs between each time point to zero time point (t0) for each sample.

	Args:
		arr_exp: 2D gene expression matrix file (genes * samples)
		lst_str2idx: a list of mapping sample, time-point, and replicate to column index.
			lst_str2idx[s] returns the list of time-points of sample s.
			lst_str2idx[s][t] returns the list of replicates of time-points t of sample s.
			lst_str2idx[s][t][r] returns the column index of replicate r of time-point t of sample s.
		lst_gene: a list of gene ids.
		valueType: types of DEG values
			'absolute_binary': 1(DEG, i.e. adj.p < 0.05) or 0(non-DEG)
			'signed_binary': 1(UP-DEG), -1(DOWN-DEG), or 0(non-DEG)
			'absolute_zstats': abs(zstats)   zstats is z where p = cdf(z) of the normal distribution cosidering UP and DOWN.
			'signed_zstats': zstats
			'absolute_FC': abs(log(fold change))
			'signed_FC': log(fold change)
		passBy: if passBy is True and the output file 'info.DEGresult' exists, then it does not run DEG detection.
		outdir: output directory where the intermediate result files are stored.
	Returns:
		arr_DEG: 2D numpy.array (genes * samples)
		lst_st2idx: a list of mapping sample and time-point to column index.
			lst_str2idx[s] returns the list of time-points of sample s .
			lst_str2idx[s][t] returns the column index of time-point t of sample s.
			Note that, DEG is a result of binary comparison of tx to t0, and t0 is excluded in arr_DEG.
			For instance, arr_DEG[g, lst_str2idx[0][0]] returns DEG value of gene g in sample s=0 of t=1
	Raises:
		pass
	"""

	if outdir == None:
		outdir='.'
	tmpfile1=outdir+'/step1.expNormalized'
	tmpfile2=outdir+'/step1.sample2group'
	tmpfile3=outdir+'/step1.compareGroup'
	tmpfile4=outdir+'/step1.DEGresult'
	if passBy == False:
		OF1=open(tmpfile1,'w')
		OF2=open(tmpfile2,'w')
		OF3=open(tmpfile3,'w')
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
		subprocess.call(["Rscript", "DEGlimma.R", tmpfile1, tmpfile2, tmpfile3, valueType, tmpfile4])
	IF = open(tmpfile4,'r')
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
	arr_DEG=np.array(lst_DEGs)
	IF.close()
	return arr_DEG, lst_st2idx

def testDEGprofile(arr_DEG, method='sampling', nTrial = 200000, seed=None):
	"""profiling DEGs between each time point to zero time point (t0) for each sample.

	Args:
		arr_exp: 2D gene expression matrix file (genes * samples)
		lst_str2idx: a list of mapping sample, time-point, and replicate to column index.
			lst_str2idx[s] returns the list of time-points of sample s.
			lst_str2idx[s][t] returns the list of replicates of time-points t of sample s.
			lst_str2idx[s][t][r] returns the column index of replicate r of time-point t of sample s.
		lst_gene: a list of gene ids.
		valueType: types of DEG values
			'absolute_binary': 1(DEG, i.e. adj.p < 0.05) or 0(non-DEG)
			'signed_binary': 1(UP-DEG), -1(DOWN-DEG), or 0(non-DEG)
			'absolute_zstats': abs(zstats)   zstats is z where p = cdf(z) of the normal distribution cosidering UP and DOWN.
			'signed_zstats': zstats
			'absolute_FC': abs(log(fold change))
			'signed_FC': log(fold change)
		passBy: if passBy is True and the output file 'info.DEGresult' exists, then it does not run DEG detection.
		outdir: output directory where the intermediate result files are stored.
	Returns:
		arr_DEG: 2D numpy.array (genes * samples)
		lst_st2idx: a list of mapping sample and time-point to column index.
			lst_str2idx[s] returns the list of time-points of sample s .
			lst_str2idx[s][t] returns the column index of time-point t of sample s.
			Note that, DEG is a result of binary comparison of tx to t0, and t0 is excluded in arr_DEG.
			For instance, arr_DEG[g, lst_str2idx[0][0]] returns DEG value of gene g in sample s=0 of t=1
	Raises:
		pass
	"""
	if method == 'sampling':
		if seed != None:
			np.random.seed(seed)
		nGene = arr_DEG.shape[0]
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
		lst_pval = np.ones(nGene)
		lst_gene_val = sorted([(g, np.sum(arr_DEG[g,:])) for g in range(nGene)],key=lambda x:x[1], reverse=True)
		lst_randomVal.sort(reverse=True)
		lst_randomVal.append(-np.inf)
		i,j=0,0
		while True:
			if i >= nGene:
				break
			if j < nTrial and lst_randomVal[j] == lst_randomVal[j+1]:
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

def multipleTestCorrection(lst_pval, method='bonferroni'):
	reject,lst_adjPval,alphacSidak,alphacBonf = statsmodels.sandbox.stats.multicomp.multipletests(lst_pval,method=method)
	return lst_adjPval

def skmeansWrapper(arr_exp, k, lst_gene=None, lst_sample=None, writeInputfile=True, passBy=False, seed=None, outdir=None):
	if outdir == None:
		outdir = '.'
	if lst_gene == None:
		lst_gene = ['g%09d'%(g) for g in range(arr_exp.shape[0])]
	if lst_sample == None:
		lst_sample = [str(i) for i in range(arr_exp.shape[1])]
	tmpfile5=outdir+'/step2.expNormalized.consensusDEG'
	tmpfile6=outdir+'/step2.clusterResult.K'+str(k)
	if passBy == False:
		if writeInputfile == True:
			OF5=open(tmpfile5,'w')
			OF5.write('ID\t'+'\t'.join(map(str,lst_sample))+'\n')
			for g in range(arr_exp.shape[0]):
				OF5.write(lst_gene[g]+'\t'+'\t'.join(map(str,arr_exp[g,:]))+'\n')
			OF5.close()
		if seed == None:
			subprocess.call(["Rscript", "run_skmeans.R", tmpfile5, str(k), 'None', tmpfile6])
		else:
			subprocess.call(["Rscript", "run_skmeans.R", tmpfile5, str(k), str(seed), tmpfile6])
	IF = open(tmpfile6,'r')
	dic_cluster2idxs={}
	dic_gene2idx={}
	for g in range(arr_exp.shape[0]):
		dic_gene2idx[lst_gene[g]]=g
	for line in IF:
		s=line.strip().split('\t')
		gene, cluster = s[0], int(s[1])
		idx=dic_gene2idx[gene]
		if cluster not in dic_cluster2idxs:
			dic_cluster2idxs[cluster]=set()
		dic_cluster2idxs[cluster].add(idx)
	IF.close()
	return dic_cluster2idxs

def getResponseTimeVector(arr_exp, lst_str2idx):
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

def run(lst_sample, lst_gene, arr_exp, DEGpcut, minClusterSize=3, step2_k='auto', platform='microarray', seed=None, outdir=None):
	########################
	#        STEP 1        #
	########################
	# Convert array to sarr
	lst_str2idx = label2struct(lst_sample)
	'''
		A lst_str2idx maps sample, time-point, replicate to column index.
			lst_str2idx[s] returns the list of time-points of sample s.
			lst_str2idx[s][t] returns the list of replicates of time-points t of sample s.
			lst_str2idx[s][t][r] returns the column index of replicate r of time-point t of sample s.
		Thus, arr_exp[g][lst_str2idx[s][t][r]] returns the the expression value of gene g in sample s, timepoint t, replicate r.
	'''
	# Normalize expression
	arr_exp = normailze(arr_exp, lst_str2idx, platform=platform)

	# Profile DEG for t0 vs. tx for each sample
	arr_DEGpre, lst_st2idx = profileDEG(arr_exp, lst_str2idx, lst_gene=lst_gene, valueType='signed_zstats', passBy=False, outdir=outdir)
	'''
		A arr_DEGpre is a numpy array of two dimensions.
		That is, arr_DEGpre[g,lst_st2idx[s][t-1]] returns whether gene g is DEG or not for timepoint t vs. t0 in sample s.
		Also valueType determines the meaning of DEG value
			'absolute_binary': 1(DEG, i.e. adj.p < 0.05) or 0(non-DEG)
			'signed_binary': 1(UP-DEG), -1(DOWN-DEG), or 0(non-DEG)
			'absolute_zstats': abs(zstats)    where zstats = cdf.inverse(p) of the normal distribution cosidering UP and DOWN.
			'signed_zstats': zstats
			'absolute_FC': abs(log(fold change))
			'signed_FC': log(fold change)
	'''

	# Merge DEG profile for each sample
	arr_DEG=np.zeros((arr_DEGpre.shape[0], len(lst_st2idx)))
	for s in range(len(lst_st2idx)):
		lst_i = lst_st2idx[s]
		arr_DEG[:,s] = np.max(abs(arr_DEGpre[:,lst_i]),axis=1)
	'''
		A arr_DEG is a numpy array of two dimensions.
		That is, arr_DEG[g,s] returns whether gene g is DEG or not in sample s.
	'''

	# Compute p-value of DEG frequencies against the randomly sampled frequency distribution.
	lst_pval = testDEGprofile(arr_DEG, method='sampling', nTrial=200000, seed=seed)

	# perform multiple test correction to produce adjusted p-values
	_a_, lst_adjPval, _b_, _c_ = statsmodels.sandbox.stats.multicomp.multipletests(lst_pval,method='fdr_bh')

	# select indices of consensus DEG genes
	lst_DEGidx = [i for i in range(len(lst_gene)) if lst_adjPval[i] < DEGpcut]

	# sort DEG genes by frequency
	lst_frequency = np.sum(arr_DEG, axis=1)
	lst_DEGidx = sorted(lst_DEGidx, key=lambda x:lst_frequency[x], reverse=True)
	lst_DEGgene = [lst_gene[i] for i in lst_DEGidx]
	arr_DEGexp = arr_exp[lst_DEGidx,:]

	# select 10% of DEG genes as pseudo reference genes
	set_pseudoRefGene = set(lst_DEGgene[0:int(0.1*len(lst_DEGgene))])

	# set the number of clusters k
	if step2_k == 'auto':
		maxK = min(500,len(arr_DEGexp))
		lst_k = list(range(50,maxK,50))
	else:
		lst_k = [int(step2_k)]

	dic_k2f1_gene2phase_phase2rt={}
	for kidx, k in enumerate(lst_k):
		########################
		#        STEP 2        #
		########################
		if kidx == 0:
			dic_cluster2idxs = skmeansWrapper(arr_DEGexp, k, lst_gene=lst_DEGgene, lst_sample=lst_sample, writeInputfile=True, seed=seed, outdir=outdir) 
		elif kidx >= 1:
			dic_cluster2idxs = skmeansWrapper(arr_DEGexp, k, lst_gene=lst_DEGgene, lst_sample=lst_sample, writeInputfile=False, seed=seed, outdir=outdir) 
		for cluster in dic_cluster2idxs.keys():
			if len(dic_cluster2idxs[cluster]) < minClusterSize:
				del dic_cluster2idxs[cluster]
	
		########################
		#        STEP 3        #
		########################
		dic_cluster2rt_tstats={}
		for i, (cluster, gidxs) in enumerate(dic_cluster2idxs.items()):
			rt, tstats = getResponseTimeVector(arr_DEGexp[sorted(gidxs),:], lst_str2idx)
			dic_cluster2rt_tstats[cluster]=[rt,tstats]

		########################
		#        STEP 4        #
		########################
		dic_phase2rt={}
		dic_gene2phase={}
		dic_cluster2phase = orderingCluster(dic_cluster2rt_tstats)
		for cluster, phase in dic_cluster2phase.items():
			for i in dic_cluster2idxs[cluster]:
				gene = lst_DEGgene[i]
				dic_gene2phase[gene]=phase
				dic_phase2rt[phase]=dic_cluster2rt_tstats[cluster][0]

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
		dic_k2f1_gene2phase_phase2rt[k]=[f1, dic_gene2phase, dic_phase2rt]
	k, (f1, dic_gene2phase, dic_phase2rt) = sorted(dic_k2f1_gene2phase_phase2rt.items(),key=lambda x:x[1][0], reverse=True)[0]
	return dic_gene2phase, dic_phase2rt


if __name__ == "__main__":
	parser=argparse.ArgumentParser(
		usage='''\
	%(prog)s [options] expfile -o out.txt
	example: %(prog)s expfile -o out.txt
	''')
	
	parser.add_argument('expfiles', nargs='+', help='Gene expression file')
	parser.add_argument('-DEGpcut', required=False, type=float, default=0.05, help='Cutoff value of adjusted p-value for testing consensus DEGs')
	parser.add_argument('-minClusterSize', required=False, type=int, default=3, help='Cutoff value of proportion for the minimum cluster size')
	parser.add_argument('-step2_k', required=False, default=200, help='The number of clusters for Step 2')
	parser.add_argument('-seed', required=False, type=int, default=None, help='Seed value to fix random process')
	parser.add_argument('-platform', required=False, choices=['microarray','RNA-seq'], default='microarray', help='The platform of measurment of gene expression data')
	parser.add_argument('-o', required=False, metavar='str', default='result', help='Output directory')
	args=parser.parse_args()
	
	if not os.path.exists(args.o):
		os.makedirs(args.o)

	lst_sample, lst_gene, arr_exp = readData(args.expfiles)

	dic_gene2phase, dic_phase2rt = run(lst_sample=lst_sample, lst_gene=lst_gene, arr_exp=arr_exp, DEGpcut=args.DEGpcut, minClusterSize=args.minClusterSize, step2_k=args.step2_k, platform=args.platform, seed=args.seed, outdir=args.o)
	OF=open(args.o+'/final.gene2phase.txt','w')
	for gene, phase in dic_gene2phase.items():
		OF.write('\t'.join(map(str,[gene, phase]))+'\n')
	OF.close()

	OF=open(args.o+'/final.phase2responseTime.txt','w')
	for phase, rt in dic_phase2rt.items():
		str_rt = ','.join(map(str,rt))
		OF.write('\t'.join(map(str,[phase, str_rt]))+'\n')
	OF.close()
