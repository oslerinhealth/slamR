**slamR**: **S**tructured **L**atent **A**ttribute **M**odels in **R**

An R package for Structured Latent Attribute Models (SLAM)

**Maintainer**: Zhenke Wu, zhenkewu@umich.edu

**References**: If you are using **slamR** for clustering multivariate binary
observations with SLAM, please cite the following paper:

|       | Citation     | Paper Link
| -------------  | -------------  | -------------  |
| SLAM - known Q    | Gu Y and Xu G (2019). Learning attribute patterns in high-dimensional structured latent attribute models. *Journal of Machine Learning Research* 20.115: 1-58   |[Link](http://jmlr.org/papers/v20/19-197.html)| 


## Table of content
- [1. Installation](#id-section1)
- [2. Overview](#id-section2)
- [2. Example](#id-section3)

<div id='id-section1'/>

Installation
--------------
```r
install.packages("devtools",repos="https://cloud.r-project.org")
devtools::install_github("zhenkewu/slamR")
```
<div id='id-section2'/>

Overview
----------
Structured latent attribute models (SLAMs) are a special family of discrete latent variable models widely used in social and biological sciences. This paper considers the problem of learning significant attribute patterns from a SLAM with potentially high-dimensional configurations of the latent attributes. We address the theoretical identifiability issue, propose a penalized likelihood method for the selection of the attribute patterns, and further establish the selection consistency in such an overfitted SLAM with a diverging number of latent patterns. The good performance of the proposed methodology is illustrated by simulation studies and two real datasets in educational assessments.

**slamR** works for 

* one-level binary responses
	-  known Q, unknown attribute set
    -  unknown Q, unknown attribute set
* two-level binary responses
	-  unknown or known Q, unknown attribute sets at both levels


<div id='id-section3'/>

Examples (two-level binary responses)
---------

- 1. restricted latent class analysis with **pre-specified # of factors** but unknown # of clusters

* Example of the Tree Structure of Observed ICD-9 Codes
![](inst/example_figure/tree.png)

* Incorporate tree structure
![](inst/example_figure/tree_D.png)

* doubly-multi-resolution approach for dealing with two-level multivariate binary data
![](inst/example_figure/model_structure.png)



