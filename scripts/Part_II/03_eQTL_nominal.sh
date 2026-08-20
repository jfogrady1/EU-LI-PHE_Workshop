# ---------------------------------------------------------------------------
# Define variables
# ---------------------------------------------------------------------------

BASE_DIR="/home/workspace/jogrady/other_projects/EU-LI-PHE_Workshop"
DATA_DIR="${BASE_DIR}/data/Part_II"
RESULTS_DIR="${BASE_DIR}/results/Part_II"
 
# input: phenotype bed (+ index) and covariates
BED_PHENOTYPE="${DATA_DIR}/ALL.expr_tmm_inv.Chr5.bed.gz"
BED_PHENOTYPE_TBI="${DATA_DIR}/ALL.expr_tmm_inv.Chr5.bed.gz.tbi"
COVARIATES_FILE="${DATA_DIR}/ALL.covariates.txt"
 
# params: genotype plink prefix
PLINK_PREFIX_PATH="${DATA_DIR}/ALL_genotypes_Chr5"
 
# output: cis-QTL results
PARQUET_FILE="${RESULTS_DIR}/ALL.cis_qtl_pairs.5.parquet"

OUTPUT_PREFIX="${RESULTS_DIR}/ALL"


# ---------------------------------------------------------------------------
# Run tensorQTL (cis, permutation pass)
# ---------------------------------------------------------------------------
echo "Running tensorQTL cis-permutation mapping"
echo "  Genotypes (plink prefix): ${PLINK_PREFIX_PATH}"
echo "  Phenotypes (bed):         ${BED_PHENOTYPE}"
echo "  Covariates:               ${COVARIATES_FILE}"
echo "  Output prefix:            ${OUTPUT_PREFIX}"
echo 
    python3 -m tensorqtl \
    "${PLINK_PREFIX_PATH}" \
    "${BED_PHENOTYPE}" \
    "${OUTPUT_PREFIX}" \
    --mode cis_nominal \
    --window 1000000 \
    --covariates "${COVARIATES_FILE}" \

echo
echo "Done. Expected output: ${PARQUET_FILE}"

OUTPUT_FILE="${RESULTS_DIR}/ALL.cis_qtl_pairs.5.txt.gz"

INPUT_SCRIPT="${BASE_DIR}/scripts/Part_II/parquet2txt.py"

echo "Converting parquet to .txt"

python3 ${INPUT_SCRIPT}  ${PARQUET_FILE} ${OUTPUT_FILE}

echo "Done. Expected output: ${OUTPUT_FILE}"