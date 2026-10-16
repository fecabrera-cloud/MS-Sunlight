This repository contains the processed data, code, tables, and manuscript figures for the accompanying study of direct and litter-mediated sunlight effects on maize-litter microbiomes.

The package begins at the processed bioinformatics-table level. Raw FASTQ processing, read trimming, host filtering, HUMAnN execution, and HPC job orchestration are intentionally outside its scope. 

Running the package reproduces:

- field and microcosm amplicon ordinations, PERMANOVA/dispersion tests, composition summaries, and genus-level differential-abundance analyses;
- field and microcosm targeted KO functional-group ordinations, univariate models, sensitivity analyses, exploratory prevalent-KO screens, and supporting plots;
- microcosm dry-mass loss, microbial biomass C, beta-glucosidase, and biomass-to-respiration ratio models;
- descriptive substrate-chemistry tables, without inferential pseudoreplication;
- field temperature and relative-humidity summaries, GEE analyses, and plots;
- Supplementary Tables S1-S3 and the complete KO curation workbook.


## Repository structure

```text
data/                    Processed input tables and sample ledger
scripts/                 Complete analysis and validation code
figures/                 Manuscript figures
tables/                  Manuscript-ready supplementary tables
environment/             Software versions and installation notes
docs/                    Analysis boundary, data dictionary, and maps
```

## Run the complete analysis

The tested environment used R 4.1.2 and Python 3.13.2. Install the packages listed in `environment/`, then run:

```bash
bash scripts/run_all.sh
```

The runner executes all analyses in dependency order and finishes with structural and numerical validation. A successful run ends with `All release checks passed.` and writes `validation/VALIDATION_REPORT.txt`.

Individual components can also be run directly from any working directory because every script resolves paths relative to its own location. See [docs/FIGURE_AND_TABLE_MAP.md](docs/FIGURE_AND_TABLE_MAP.md).

## Treatment terminology

- `S+` / `SusS+`: full sunlight exposure or litter preconditioned under full sunlight.
- `UV-` / `SusUV-`: reduced UV, with visible wavelengths retained.
- `S-` / `SusS-`: reduced UV and blue-green light (strongest spectral attenuation).

The `Sus` prefix is used for microcosm substrate photo-history. Analyses of the two biological experiments remain separate.

## Functional-metagenomic scope

The primary functional analysis is targeted and hypothesis driven. It uses a frozen mapping of 332 KEGG Orthology identifiers to 13 mutually exclusive functional groups, drawn from a treatment-blind universe of 7,931 detected KOs. These groups represent approximately 5% of total KO relative abundance and are not a complete description of community metabolism. Each KO is assigned to one primary group to prevent double counting; the full evidence and exclusion audit is included under `data/functional_metagenomics/` and Table S2.

HUMAnN profiles are unstratified, so they do not identify the organisms carrying a function. Amplicon-responsive taxa therefore cannot be directly assigned responsibility for KO-group patterns.

## Licensing

Code is released under the MIT License (`LICENSE-CODE`). Data, tables, documentation, and figures are released under Creative Commons Attribution 4.0 International (`LICENSE-DATA`). Third-party databases retain their own licenses and are not redistributed here.
