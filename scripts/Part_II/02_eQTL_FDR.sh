BASE_DIR="/home/workspace/jogrady/other_projects/EU-LI-PHE_Workshop"
DATA_DIR="${BASE_DIR}/data/Part_II"
RESULTS_DIR="${BASE_DIR}/results/Part_II"

PERM_RESULTS="${RESULTS_DIR}/ALL.cis_qtl.txt.gz"
INPUT_SCRIPT="${BASE_DIR}/scripts/Part_II/eGene_FDR_calculation.R"

CORRECTED_RESULTS="${RESULTS_DIR}/ALL.cis_qtl_fdr0.05.txt"

# Here, we are calling an R script inside a bash script
# Each argument is passed to the R script as a command line argument
# argument 1 [${INPUT_FILE}] on the command line corresponds to ARGS[1] in the R script, and so on
# argument 2 is output file (not gzipped)
# argument 3 is the FDR threshold (0.05)
echo 'Calculating global FDR for cis-eQTL permutation results'
Rscript "${INPUT_SCRIPT}" "${PERM_RESULTS}" "${CORRECTED_RESULTS}" 0.05