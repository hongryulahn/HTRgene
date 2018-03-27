# HeteroTimeDEG
This program was designed for integrating analysis of multiple heterogeneous time-series gene expression data to identify response genes


## Installation

To download all the examples, simply clone this repository:
```
git clone https://github.com/hongryulahn/HeteroTimeDEG
```

## Dependancy
To run them, you also need the latest version of the dependency library below.

### python 
- **numpy, scipy** 
- **[statsmodels](https://www.statsmodels.org/stable/install.html)**

### R
- **[skmeans](https://cran.r-project.org/web/packages/skmeans/index.html)**
- **[limma](https://bioconductor.org/packages/release/bioc/html/limma.html)**

## Input file format
Multiple time-series gene expression data has to be stored into a single tab-dilimenated matrix file (genes * samples).

### Sample labels
The first row of the file includes labels of samples in the format of i_j_k.
- i : ith time-series sample
- j : jth time-point in the sample i
- k : kth replicate of the jth time-point in the sample i

### Gene IDs
The first column of the file includes gene IDs.

## Run
```
python heteroTimeDEG.py 

#### 5 - Data Management
- **Build an image dataset** ([notebook](https://github.com/aymericdamien/TensorFlow-Examples/blob/master/notebooks/5_DataManagement/build_an_image_dataset.ipynb)) ([code](https://github.com/aymericdamien/TensorFlow-Examples/blob/master/examples/5_DataManagement/build_an_image_dataset.py)). Build your own images dataset with TensorFlow data queues, from image folders or a dataset file.
- **TensorFlow Dataset API** ([notebook](https://github.com/aymericdamien/TensorFlow-Examples/blob/master/notebooks/5_DataManagement/tensorflow_dataset_api.ipynb)) ([code](https://github.com/aymericdamien/TensorFlow-Examples/blob/master/examples/5_DataManagement/tensorflow_dataset_api.py)). Introducing TensorFlow Dataset API for optimizing the input data pipeline.
