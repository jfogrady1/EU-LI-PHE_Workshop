#!/usr/bin/env bash
#!/usr/bin/env bash

set -euo pipefail

# ==============================================================================
# Directories
# ==============================================================================

BASE_DIR="/home/workspace/jogrady/other_projects/EU-LI-PHE_Workshop"
DATA_DIR="${BASE_DIR}/data/Part_II"
RESULTS_DIR="${BASE_DIR}/results/Part_II"

mkdir -p "${RESULTS_DIR}"

# ==============================================================================
# Input files
# ==============================================================================

# Phenotype BED
BED_PHENOTYPE="${DATA_DIR}/ALL.expr_tmm_inv.Chr5.bed.gz"

# Covariates
COVARIATES_FILE="${DATA_DIR}/ALL.covariates.txt"

# Interaction term (sample IDs in first column)
INTERACTION_FILE="${DATA_DIR}/ALL_TB_interaction_input.txt"

# PLINK genotype prefix
PLINK_PREFIX_PATH="${DATA_DIR}/ALL_genotypes_Chr5"

# ==============================================================================
# Output
# ==============================================================================

OUTPUT_PREFIX="${RESULTS_DIR}/ALL_TB_interaction"

# ==============================================================================
# Run interaction eQTL mapping
# ==============================================================================


# Run TensorQTL interaction eQTL mapping
python3 -m tensorqtl \
    "${PLINK_PREFIX_PATH}" \
    "${BED_PHENOTYPE}" \
    "${OUTPUT_PREFIX}" \
    --covariates "${COVARIATES_FILE}" \
    --window 1000000 \
    --seed 1864 \
    --interaction "${INTERACTION_FILE}" \
    --maf_threshold_interaction 0.05 \
    --mode cis_nominal

echo "Interaction eQTL analysis complete."
echo "Results written to: ${OUTPUT_PREFIX}.cis_qtl_top_assoc.txt.gz"