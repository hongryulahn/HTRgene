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

Thus, the label of a column that belongs to time-series sample 0, time-point 0, replicate 2 will be 0_0_2.

### Gene IDs
The first column of the file includes gene IDs.

see the example file



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
