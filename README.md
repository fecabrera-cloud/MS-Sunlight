# Sunlight and maize-litter microbiomes

Processed scientific data, analysis code, tables, and figures for the study of direct and litter-mediated sunlight effects on maize-litter microbiomes.

The package begins at the processed bioinformatics-table level. It reproduces field and microcosm amplicon analyses, targeted functional profiles, microcosm functional-group and measured-function factorial models, substrate chemistry summaries, and field microclimate analyses. The genus-level DESeq2 analysis retains its existing specification; see [statistical methods](docs/ANALYSIS_METHODS.md).

## Structure

- `data/`: processed inputs, frozen KO curation, ordination reference inputs, numerical validation targets, and derived field abundance tables.
- `scripts/`: analysis, rendering, table assembly, and validation.
- `figures/final/`: Figures 1–6 and S1–S11, with a derived supplementary figure bundle.
- `tables/`: supplementary Tables S1–S2 and the complete KO mapping workbook.
- `environment/`: installation guidance and tested versions.
- `docs/`: methods, data dictionary, figure map, and asset checksums.

## Reproduce

Install the dependencies in `environment/`, then run:

```bash
bash scripts/run_all.sh
```

The workflow writes `results/reproduced/` and `validation/`, which are not tracked. It ends with numerical and structural release checks. Final PDF/PPTX figure assets are not overwritten. Individual scripts resolve paths relative to the repository and can be launched from any working directory.

Microcosm community tests use fixed inoculum, substrate, and interaction terms with within-block permutations. Microcosm functional-group and measured-function models retain random block. The prevalence rule is `max(2, ceiling(0.05 * N))` before rarefaction. The targeted KO analysis uses 332 mutually exclusive assignments across thirteen predefined groups, approximately 5% of total KO abundance, and does not describe the complete metagenome. HUMAnN profiles are unstratified and do not identify the organisms carrying each function.

## Treatment terminology

`S+`/`SusS+` denote full sunlight, `UV-`/`SusUV-` reduced UV, and `S-`/`SusS-` reduced UV and blue-green. `Inoc` prefixes identify inoculum origin. The `Sus` prefix identifies litter used as microcosm substrate. Figures use treatment names; retained genus-volcano contrasts also use the internal codes defined above.

## Licensing

Code: MIT (`LICENSE-CODE`). Data, documentation, tables, and figures: CC BY 4.0 (`LICENSE-DATA`). Third-party databases retain their own licenses and are not redistributed here.
