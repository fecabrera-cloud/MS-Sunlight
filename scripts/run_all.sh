#!/usr/bin/env bash
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"
mkdir -p results/reproduced validation tables
FIGURE_HASHES_BEFORE="$(mktemp)"
FIGURE_HASHES_AFTER="$(mktemp)"
trap 'rm -f "$FIGURE_HASHES_BEFORE" "$FIGURE_HASHES_AFTER"' EXIT
find figures/final -type f \( -name '*.pdf' -o -name '*.pptx' \) -print0 | sort -z | xargs -0 sha256sum > "$FIGURE_HASHES_BEFORE"
echo '[1/8] Amplicon preprocessing and community models'
Rscript scripts/run_amplicon_audit.R
echo '[2/8] Genus-level differential abundance'
python scripts/run_amplicon_deseq.py
echo '[3/8] Functional profiles and functional groups'
Rscript scripts/run_functional_metagenomics_audit.R
echo '[4/8] Microcosm functional measurements'
Rscript scripts/run_functional_measurements_audit.R
echo '[5/8] Field and microcosm genus comparison'
Rscript scripts/plot_genus_comparison.R
echo '[6/8] Substrate chemistry summaries'
python scripts/run_chemistry_analysis.py
echo '[7/8] Field microclimate'
python scripts/run_microclimate_analysis.py
echo '[8/8] Supplementary tables, figure bundle, and validation'
python scripts/build_supplementary_tables.py
python scripts/build_figure_bundle.py
find figures/final -type f \( -name '*.pdf' -o -name '*.pptx' \) -print0 | sort -z | xargs -0 sha256sum > "$FIGURE_HASHES_AFTER"
cmp "$FIGURE_HASHES_BEFORE" "$FIGURE_HASHES_AFTER"
python scripts/validate_release.py
