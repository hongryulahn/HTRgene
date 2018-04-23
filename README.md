# HTRgene
This program was designed for integrating analysis of multiple heterogeneous time-series gene expression data to identify response genes



## Installation

To download all the examples, simply clone this repository:
```
git clone https://github.com/hongryulahn/HTRgene
```


## Dependancy
To run them, you also need the latest version of the dependency library below.

### python 
- **[numpy & scipy](https://www.scipy.org/install.html)**
- **[statsmodels](https://www.statsmodels.org/stable/install.html)**

### R
- **[skmeans](https://cran.r-project.org/web/packages/skmeans/index.html)**
- **[limma](https://bioconductor.org/packages/release/bioc/html/limma.html)**



## Input file format
Multiple time-series gene expression data has to be stored into a single tab-dilimenated matrix file (genes * samples).

### Sample labels
The first line of the file includes labels of columns in the format of i_j_k.
- i : time-series sample
- j : time-point
- k : replicate

For instance, the label of a column 0_0_2 represents time-series sample 0, time-point 0, replicate 2. 

### Gene IDs
The first column of the file includes gene IDs.

### Example
Below is an example of two time series samples, of four time points, with two replicates, for four genes.
```
ID	0_0_0	0_0_1	0_1_0	0_1_1	0_2_0	0_2_1	0_3_0	0_3_1	1_0_0	1_0_1	1_1_0	1_1_1	1_2_0	1_2_1	1_3_0	1_3_1
AT1G01010	5.35146	5.215	4.54896	4.69377	5.27511	4.53756	4.50158	4.57903	6.94247	6.4214	6.74755	6.9412	6.94886	7.06849	7.70403	7.27114
AT1G01030	5.00595	5.17827	4.96035	5.21129	5.38342	5.08703	5.01766	5.3284	4.19308	4.16198	4.16957	4.32094	4.25848	4.4901	4.31466	4.2253
AT1G01040	7.54401	7.41923	7.19867	7.16885	7.23367	7.43978	7.49953	7.57841	7.55364	7.54586	7.55467	7.37159	7.44254	7.15361	7.51412	7.71205
AT1G01050	9.52625	9.71189	9.62682	9.74527	9.85147	9.70006	9.56717	9.34816	9.76466	10.0955	10.1058	9.92132	10.0833	9.97274	9.83797	10.0723
```



## Run HTRgene
```
python HTRgene.py example/exp.cold.28sample.txt -o out.cold
python HTRgene.py example/exp.heat.24sample.txt -o out.heat
```

## Plot heatmap
```
python Heatmap.py out.cold/step1.DEGresult out.cold/final.gene2phase.txt
python Heatmap.py out.heat/step1.DEGresult out.heat/final.gene2phase.txt
```
