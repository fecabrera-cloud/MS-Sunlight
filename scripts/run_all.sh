#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

mkdir -p results/reproduced validation tables

# The manuscript-ready figures are author-controlled assets. Analyses write only
# to results/reproduced; this checksum snapshot enforces that boundary.
FIGURE_HASHES_BEFORE="$(mktemp)"
trap 'rm -f "$FIGURE_HASHES_BEFORE"' EXIT
find figures/final -type f -print0 | sort -z | xargs -0 sha256sum > "$FIGURE_HASHES_BEFORE"

echo "[1/8] Amplicon community analysis"
Rscript scripts/run_amplicon_audit.R

echo "[2/8] Amplicon genus-level differential abundance"
python scripts/run_amplicon_deseq.py

echo "[3/8] Targeted functional metagenomics"
Rscript scripts/run_functional_metagenomics_audit.R

echo "[4/8] Functional-group sensitivity analysis"
Rscript scripts/run_functional_group_sensitivity.R

echo "[5/8] Microcosm functional measurements"
Rscript scripts/run_functional_measurements_audit.R

echo "[6/8] Descriptive substrate chemistry"
python scripts/run_chemistry_analysis.py

echo "[7/8] Field microclimate"
python scripts/run_microclimate_analysis.py

echo "[8/8] Supplementary tables and release validation"
python scripts/build_supplementary_tables.py

FIGURE_HASHES_AFTER="$(mktemp)"
find figures/final -type f -print0 | sort -z | xargs -0 sha256sum > "$FIGURE_HASHES_AFTER"
if ! cmp -s "$FIGURE_HASHES_BEFORE" "$FIGURE_HASHES_AFTER"; then
  echo "ERROR: analysis modified an author-controlled file in figures/final" >&2
  rm -f "$FIGURE_HASHES_AFTER"
  exit 1
fi
rm -f "$FIGURE_HASHES_AFTER"

python scripts/validate_release.py
echo "Complete. See validation/VALIDATION_REPORT.txt."
