# EU-LI-PHE Functional Genomics Workshop
Repository for functional genomics practical component of EU-LI-PHE COST Action Summer School on Epidemiology, Genetics and Modelling for Infectious Disease Control ​

## Tutors

**Lead Tutor:** Dr. John F. O'Grady (ETH Zurich; jogrady@ethz.ch)  

**Co-Tutor:** Prof. David MacHugh (University College Dublin; david.machugh@ucd.ie)

## Learning Outcomes

The aim of this workshop is to gain practical experience in performing multi-omic data analysis to understanding bovine tuberculosis (bTB) disease. The specific learning outcomes are as follows:

- Understand how population structure can be inferred and visualized from high-resolution genomic data and appreciate why it is important to account for this variation.

- Learn how to conduct a pairwise differential expression analysis in R, accounting for potential confounders.

- Learn how to appropriately conduct a gene-set enrichment/overrepresentation analysis of differentially expressed genes (DEG).

- Understand the statistical methodologies underpinning cis-eQTL mapping and learn how to conduct an integrative molecular QTL mapping analysis.

- Become familiar with visualising outputs of an integrative genomics analysis in R

- Obtain an appreciation for how one appropriately adjusts for multiple hypothesis tests to mitigate against type I errors (false positives).


## Structure of Workshop
This workshop has been divided up into two parts, each with four exercises:

### Part I - Population genetic structure and Differential Expression Analysis (DEA). 

1. Visualization and interpretation of interferon gamma release assay (IGRA) test results.

2. Visualization and interpretation of population genetic structure in livestock.

3. Differential expression analysis (DEA)

4. Overepresentation or gene-set enrichement analysis analysis (GSEA)

### Part II - _cis_-Expression quantitative trait loci (_cis_-eQTL) mapping. 

1. Expression quantitative trait loci (eQTL) mapping (in _cis_)

2. Implementation of appropriate multiple testing correction.

3. Visualiation of _cis_-eQTL properties.

4. Interaction ieQTL mapping. 


## The data
To understand these methodologies, we will focus on bovine tuberculosis (bTB), an endemic disease of cattle caused by infection with _Mycobacterium bovis_. We will be using an updated and smaller variant of the dataset published from [O'Grady et al., (2025)](https://doi.org/10.1038/s42003-025-07846-x)

Specficially, the dataset consists of the following:
1. Transcriptomics data from _n_ = 63 control non-infected cattle (bTB-) and _n_ = 60 _M. bovis_-infected cattle.

2. Imputed and filtered WGS data for all _n_ = 123 on chromosome 5 (BTA5)



## Installation instructions for Part I

All analysis for Part I will be conducted in R/Rstudio. This makes it relatively straight forward to get set up as you should have some packages installed from Day 1 and Day 2. However, you will need to install the following packages to complete the workshop.

### CRAN packages
You can install the following cran packages with this command in R
- `install.packages(c('tidyverse', 'data.table', 'ggplot2', 'patchwork', 'dplyr', 'gprofiler2', 'ggrepel', 'viridis', 'RcolorBrewer', 'BiocManager'))`

### Bioconductor packages
- `BiocManager::install(c("DESeq2", 'apeglm'))`

- In the scripts, these packages are called when requried with the `library()` function.

## Installation instructions for Part II

### Windows Subsystem for Linux intsallation for Windows users

Approximately 85% of this workshop will be conducted through R/R-Studio, in keeping with the theme of the Summer School. However, there are some tasks that require the use of the command line interface (CLI) or the terminal.

*Those with Apple/Mac devices. or with WSL already installed can skip this step and procede to 'Python3.11, pip3.11, python3.11-venv, and tensorQTL installation' as the terminal is readily available on these machines*.

For Windows users, you will need to install a linux interface so that you can use the terminal. This can be done using this link: https://learn.microsoft.com/en-us/windows/wsl/install 

 - Open powershell in administrator mode

 - Type the following command: `wsl --install`

 - Set a name (just type it) for your account when it asks: 'Create a default Unix user account:'

 - Set a password by typing it (you won't be able to see the outputs of this so make sure you remember it)
 
 - Type `y` or `n` for data collection question

 - Please execute this command: `sudo apt-get update` (You will need your password that you set above); This installs some packages/modules/functions that we need to install our package of interest later

 - Please also run the following: `sudo apt install build-essential`

 - Please also run the following: `sudo apt install python3-dev zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev`
 
The installation is very straightfoward. Once installed you have a linux terminal on your windows machine. Here you can run and install any linux program. For small data sets this is almost certainly more than enough.

Once the installation is completed, you can close powershell. To relaunch WSL you can
1. Open powershell in regular mode and type: `wsl.exe`

2. Search for 'ubuntu' in your machine search bar and you should see an application that you can open

Getting experience in the terminal is very important in Genomics Data Science, and in any data science discipline. This is becoming even more pertinent because massive amounts of data can be readily generated generated. In many cases remote high-performance computing (HPC)/cloud computing (e.g., Amazon Web Services AWS) servers are required to store this data. These are often linux/Ubuntu based and require some knowledge of the terminal.

Being trained in the terminal is beyond the scope of this course (We will only use it to execute bash scripts in Part_II) but there are some good resources here:
 - https://www.youtube.com/watch?v=v392lEyM29A  (Nice video tutorial on navigating the terminal)

 - https://www.geeksforgeeks.org/linux-unix/linux-commands-cheat-sheet/ (Nice cheat sheet)

Don't worry if you cannot install WSL (e.g., due to administrative privilege issues); all output files will be made available so you can continue with the practical as we go through it together!


### Python3, pip3, python3-venv, and TensorQTL installation

Once WSL installation is completed (or if you didn't need to install it), we need to install python3.

You may have python3 already installed on your Windows machine but because linux applications use a different compiled version of python, we may need to install it within this (WSL) application

- Please check that python3.11 is installed by typing: `python3.11` (or default `python`)

- If it executes, you are good to go! If it says something like 'python3.11 not found', please install python3.11 with the following command: `sudo apt install python3.11 python3.11-dev` (You may need your password for this that you set previously)


Now we need to install pip3, which is required to install tensorQTL that we will use for the molecular QTL mapping.
- Please check that you have pip3 installed by typing: `pip3 --version`

- If it executes, you are good to go! If it says something like 'pip3 not found or some other version', please install it with the following command: `sudo apt install python3.11-pip` (You may need your password for this that you set previously). Type 'Y' when prompted


Now we need to install `python3.11-venv` - this will create a virtual environment that will host all the dependancies associated with tensorqtl (essentially to make it work ok)

- Please type the following command: `sudo apt install python3.11-venv` (You may need your password) Press Y when prompted

- Now we need to create the environment and we will do so by typing the following command: `python3.11 -m venv tensorqtl` (this might take two or three minues)

- Now we need to activate the environment and we will do so with the following command: `source tensorqtl/bin/activate` (you should see 'tensorqtl' appear at the start of your prompt - you are now in this environment)

Now we need to install tensorqtl: (https://github.com/broadinstitute/tensorqtl).

- Please type the following command: `pip3 install tensorqtl` (This will take about 20 mins)

- You will next need to install pandas-plink: `pip3 install pandas-plink`

### R on linux and r2py python packages
The last thing we need for tensorQTL to work is R (as it internally also calls some R packages)
You may have R already installed on your Windows machine but because linux applications use a different compiled version of python, we may need to install it within this (WSL) application

For TensorQTL, most of you will not have this installed. 

- You can install it using the following command: 'sudo apt install r-base r-base-dev'

- You can install the r2py python package using the following command: pip install rpy2

### Specific R package for tensorqtl on WSL
On the terminal execute the following. Remember, you probably have BiocManager installed on your windows machine but you need to also install it for this compiled version of R for linux.
- Type the command: `R`
- Within R, type and execute the command `install.packages(BiocManager)`
- Type `yes` when prompted and `yes` again
- Type the following command `BiocManager::install(c('qvalue'))`
- Exit out of R with the command `ctrl Z`.

Verify that tensorqtl is installed with the following: `python3 -m tensorqtl`

*You should now be ready to do the practical*
