# Scripts for Part II of the work shop

1. `01_eQTL_permute.sh`: This script performs permutation pass in `tensorqtl` and identifies top variants per gene


2. `02_eQTL_FDR.sh`: This script calls an R script internally (see below) and corrects for multiple genes being tested gene/chromosome wide using the Benjamini-Hochberg procedure by default.
    - `eGene_FDR_calculation.R`: This scripts implements alternative multiple testing corrections chromosome/genome-wide and is called within the main bash script.
    - `eQTL_properties.R` : This script will be used to visualize and calculate the properties of _cis_-eQTLs

3. `03_eQTL_nominal.sh`: This script performs a nominal pass (i.e., associates all variants in cis using standard linear regression). Note: This file will be very long as all statistics for all variants will be output - some will have mutliple statistics if they lie in the _cis_ window of multiple genes. The ouput of this command is a `.parquet` file (alternative to .gz) file but we will convert it to .gz format.
    - `parquet2txt.py` - This is a python script that converts the parquet file to a compressed txt file (`.txt.gz`), which makes it easier to read into R etc.


4. `04_ieQTL_mapping.sh`: This script performs the interaction eQTL mapping using tensorQTL fitting an interaction term between genotype and bTB status. Same inputs as the permutation procedure but with an additional `--interaction` flag and file, respectively.
    - `ieQTL_visualization.R`: This script will be used to plot the ieQTL we identify on BTA5

    ## Important points

    1. Don't forget to activate your source environment for tensorQTL before executing the `.sh` scripts. e.g., `source ~/tensorqtl/bin/activate`.
    2. Don't forget to change your `$BASE_DIR` variable in the script to what ever location you have downloaded the repo into. e.g., I downloaded mine to my Downloads directory. `/mnt/c/Users/jogrady/Downloads/EU-LI-PHE_Workshop-main/EU-LI-PHE_Workshop-main`. And I will set that as my `$BASE_DIR`.